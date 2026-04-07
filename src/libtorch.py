# -*- coding: utf-8 -*-
"""libtorch.py — SymPy-to-Torch conversion and numerical integration helpers.

Part of libphysics.

Two-stage design for performance:
    1. libtorch.torchify()              — expensive SymPy work done ONCE
                                         (change-of-vars + ``lambdify``)
    2. TorchExpr.torchquad_integrate()  — cheap torch-only integration,
                                         called MANY times

USAGE
=====
    from libphysics.libtorch import libtorch

    # Step 1: build once (slow: SymPy subs + lambdify, ~200–500 ms)
    lt = libtorch()
    texpr = lt.torchify(
        sp.Integral(integrand, (y_A, -sp.oo, sp.oo), (y_B, -sp.oo, sp.oo))
    )

    # Step 2: integrate many times (fast: pure torch, ~3–8 ms each)
    re, im = texpr.torchquad_integrate(N=121)

    # Batched integration over a parameter grid using Simpson rule
    re_b, im_b = texpr.torch_integrate_batched_simpson(
        params_values=params_grid,
        N=21,
        chunk_size_params=256,
        chunk_size_points=10_000,
    )
"""
from dataclasses import dataclass
import sys
from typing import List, Callable, Any

import numpy as np
import torch
from scipy import special as scipy_special
from sympy import (
    Eq,
    Function,
    Integral,
    Integer,
    Symbol,
    atan,
    atanh,
    cos,
    cosh,
    diff,
    lambdify,
    oo,
    pi,
    sinh,
    tan,
    tanh,
)

from loguru import logger

def setup_logging(enable: bool = False, level: str = "INFO"):
    """
    By default, we remove all handlers (silence).
    If enable=True, we add a standard output handler.
    ------
    Example:
        # Disable logging (default)
        setup_logging()
        # Enable logging to stderr at INFO level
        setup_logging(enable=True, level="INFO")
        # Enable logging to stderr at DEBUG level
        setup_logging(enable=True, level="DEBUG")
    """
    logger.remove() 
    if enable:
        logger.add(sys.stderr, level=level)

setup_logging(enable=False) 


@dataclass
class TorchExpr:
    """Pre-built torch function together with its finite integration domain.

    This is the object returned by ``libtorch.torchify`` when integration
    limits are supplied (directly or via a SymPy ``Integral``).  It stores
    both the numerical callable and all metadata needed by the integration
    helpers in this module.

    Attributes
    ----------
    func : Callable
        The torch-compatible function produced by ``lambdify``.
        Its signature is::

            func(new_var_0, ..., new_var_{dim-1}, param_0, ..., param_{n_params-1})

        where the ``new_var_*`` correspond to the transformed integration
        variables after any change-of-variables and the ``param_*`` are any
        additional symbolic parameters.
    domain : list[list[float]]
        Finite box for the integration domain.  Each entry is ``[a, b]``
        for one dimension after change-of-variables.  This is what the
        integrators use as their integration range.
    dim : int
        Number of integration variables (i.e. the dimensionality of the
        domain).
    n_params : int
        Number of additional symbolic parameters expected by ``func``.
    sympy_expr : Any
        The SymPy expression *after* all symbolic substitutions and
        Jacobians have been applied.  Kept for debugging / inspection.
    variables : list[sympy.Symbol] | None
        Integration variables used by the numerical object after any
        change-of-variables has been applied.  For infinite limits these
        are the transformed symbols such as ``t_x`` rather than the
        original symbolic variable ``x``.
    """
    func: Callable                # f(new_var_0, …, new_var_n, param_0, …, param_m)
    domain: List[List[float]]     # finite box for torchquad
    dim: int                      # number of integration variables
    n_params: int                 # number of extra parameters
    sympy_expr: Any               # sympy expression reduced if
    variables: List[Symbol] = None       # sympy variables after change of variables

    @staticmethod
    def _rule_simpson(a: float, b: float, N: int, *, device=None, dtype=None):
        if N < 3 or (N % 2) == 0:
            raise ValueError(f"Simpson rule requires odd N>=3; got N={N}")
        coords = torch.linspace(a, b, N, device=device, dtype=dtype)
        dx = (b - a) / (N - 1)
        w = torch.ones(N, device=device, dtype=dtype)
        w[1:-1:2] = 4
        w[2:-1:2] = 2
        weights = w * (dx / 3.0)
        return coords, weights

    @staticmethod
    def _rule_trapezoid(a: float, b: float, N: int, *, device=None, dtype=None):
        if N < 2:
            raise ValueError(f"Trapezoid rule requires N>=2; got N={N}")
        coords = torch.linspace(a, b, N, device=device, dtype=dtype)
        dx = (b - a) / (N - 1)
        w = torch.ones(N, device=device, dtype=dtype)
        w[0] = 0.5
        w[-1] = 0.5
        weights = w * dx
        return coords, weights

    @staticmethod
    def _rule_gauss_legendre(a: float, b: float, N: int, *, device=None, dtype=None):
        if N < 1:
            raise ValueError(f"Gauss-Legendre rule requires N>=1; got N={N}")
        nodes, weights = np.polynomial.legendre.leggauss(N)
        nodes_t = torch.as_tensor(nodes, device=device, dtype=dtype)
        weights_t = torch.as_tensor(weights, device=device, dtype=dtype)
        # Affine map from [-1, 1] to [a, b]
        coords = 0.5 * (b - a) * nodes_t + 0.5 * (b + a)
        weights = 0.5 * (b - a) * weights_t
        return coords, weights

    @staticmethod
    def _resolve_quadrature_rule(method):
        if callable(method):
            return method

        if method is None:
            method = "simpson"

        method_name = str(method).strip().lower()
        if method_name == "simpson":
            return TorchExpr._rule_simpson
        if method_name == "trapezoid":
            return TorchExpr._rule_trapezoid
        if method_name in {"gauss-legendre", "gauss_legendre", "legendre"}:
            return TorchExpr._rule_gauss_legendre

        raise ValueError(
            "Unknown integration method. Expected one of {'simpson', 'trapezoid', 'gauss-legendre'} or a callable. "
            f"Got: {method}"
        )

    @staticmethod
    def _is_scalar_like(x) -> bool:
        """Return ``True`` for Python scalars or 0-dim tensors."""
        if torch.is_tensor(x):
            return x.ndim == 0
        return isinstance(x, (int, float, complex, np.number))

    @staticmethod
    def _normalize_params_values(params_values, n_params: int, *, device, dtype):
        """Normalize parameter inputs into a 2D tensor and batch shape.

        Accepted input formats are:

        - ``None`` when ``n_params == 0``,
        - a tensor of shape ``(..., n_params)``, or
        - a list/tuple of length ``n_params`` whose elements are scalars or
          tensors with the same shape.
        """
        if params_values is None:
            if n_params != 0:
                raise ValueError(f"Expected {n_params} params; got None")
            return torch.empty((1, 0), device=device, dtype=dtype), tuple()

        if torch.is_tensor(params_values):
            params_tensor = params_values
            if params_tensor.shape[-1] != n_params:
                raise ValueError(
                    f"params_values last dimension must be n_params={n_params}; got shape={tuple(params_tensor.shape)}"
                )
            batch_shape = tuple(params_tensor.shape[:-1])
            batch_size = int(np.prod(batch_shape)) if batch_shape else 1
            params_flat = params_tensor.reshape(batch_size, n_params).to(device=device, dtype=dtype)
            return params_flat, batch_shape

        if not isinstance(params_values, (list, tuple)):
            params_values = [params_values]

        if len(params_values) != n_params:
            raise ValueError(f"Expected {n_params} params; got {len(params_values)}")

        batch_shape = None
        for param_value in params_values:
            if torch.is_tensor(param_value) and param_value.ndim > 0:
                batch_shape = tuple(param_value.shape)
                break
            if not TorchExpr._is_scalar_like(param_value):
                param_tensor = torch.as_tensor(param_value)
                if param_tensor.ndim > 0:
                    batch_shape = tuple(param_tensor.shape)
                    break
        if batch_shape is None:
            batch_shape = tuple()

        params_list = []
        for param_value in params_values:
            param_tensor = torch.as_tensor(param_value, device=device)
            if param_tensor.ndim == 0:
                if batch_shape:
                    param_tensor = param_tensor.expand(batch_shape)
            else:
                if tuple(param_tensor.shape) != batch_shape:
                    raise ValueError(
                        f"All parameter tensors must have the same shape; got {tuple(param_tensor.shape)} vs {batch_shape}"
                    )
            params_list.append(param_tensor)

        params_tensor = torch.stack(params_list, dim=-1).to(dtype=dtype)
        batch_size = int(np.prod(batch_shape)) if batch_shape else 1
        params_flat = params_tensor.reshape(batch_size, n_params)
        return params_flat, batch_shape
    
    # ------------------------------------------------------------------
    # Convenience methods: numeric integration on this TorchExpr
    # ------------------------------------------------------------------

    def torchquad_integrate(self, params_values=None, method=None, N: int = 21):
        """Integrate this ``TorchExpr`` using the torchquad-based integrator.

        This is the simplest integration entry point.  It is intended for
        low-dimensional problems where a single torchquad integral is
        sufficient.

        Parameters
        ----------
        params_values : list | tuple | torch.Tensor | None, optional
            Numerical values for the symbolic parameters.  The accepted
            formats match those of ``torchquad_integrate``:

            - ``None`` if there are no parameters,
            - a Python sequence (list/tuple) of numbers/tensors, or
            - a tensor whose last dimension is ``n_params``.
        method : torchquad integrator instance or None, optional
            If ``None``, a default ``Simpson`` integrator from torchquad is
            constructed internally.
        N : int, optional
            Resolution parameter passed through to torchquad.  Larger values
            increase accuracy but also the number of function evaluations.

        Returns
        -------
        re, im : torch.Tensor
            Real and imaginary parts of the integral.  If the integrand is
            real-valued, ``im`` will be a zero tensor.

        USAGE
        =====
        1. Simple 1D integral without parameters::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Integral, exp, oo

            x = symbols("x", real=True)
            integral_expr = Integral(exp(-x**2), (x, -oo, oo))

            lt = libtorch()
            texpr = lt.torchify(integral_expr)
            re, im = texpr.torchquad_integrate(N=121)
            ````

        2. Integral with parameters (e.g. Fourier transform)::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Integral, exp, I, oo
            import torch

            x, k = symbols("x k", real=True)
            integral_expr = Integral(exp(-x**2) * exp(I * k * x), (x, -oo, oo))

            lt = libtorch()
            texpr = lt.torchify(integral_expr)

            # Evaluate at k = 0.5
            re, im = texpr.torchquad_integrate(params_values=[0.5], N=121)
            ````
        """
        from torchquad import Simpson

        if method is None:
            method = Simpson()

        param_vals = list(params_values) if params_values else []

        def integrand(domain_points):
            param_vals_local = [torch.as_tensor(param_value, device=domain_points.device) for param_value in param_vals]
            values = self.func(*[domain_points[:, i] for i in range(domain_points.shape[1])], *param_vals_local)
            return torch.as_tensor(values, device=domain_points.device)

        result = method.integrate(
            integrand,
            dim=self.dim,
            N=N,
            integration_domain=self.domain,
        )
        return (
            result.real if torch.is_complex(result) else result,
            result.imag if torch.is_complex(result) else torch.zeros_like(result),
        )

    def torch_integrate_batched(
        self,
        *,
        params_values = None,
        method: str | Callable = "simpson",
        N: int = 121,
        chunk_size_params: int = 256,
        chunk_size_points: int | None = None,
        device=None,
        dtype=None,
    ):
        """Batched tensor-product quadrature integration for this ``TorchExpr``.

        This method is intended for situations where you want to evaluate
        the same integral for many different parameter values (for example,
        on a 2D or 3D mesh).

        Parameters
        ----------
        params_values : tensor | list | tuple
            Numerical parameters.  See ``torch_integrate_batched``
            for a full description of the accepted shapes and types.
        method : str or callable, optional
            Quadrature rule name (``"simpson"``, ``"trapezoid"``,
            ``"gauss-legendre"``) or a callable returning
            ``(coords_1d, weights_1d)`` for ``(a, b, N)``.
        N : int, optional
            Odd number of Simpson points per dimension.
        chunk_size_params : int, optional
            Number of parameter points processed in one chunk.  Reducing
            this lowers peak memory usage at the cost of more Python loops.
        chunk_size_points : int or None, optional
            Maximum number of grid points processed per chunk.  ``None``
            means "all at once".
        device : torch.device or None, optional
            Device used for all internal tensors.  If ``None``, CUDA is
            used when available, otherwise CPU.
        dtype : torch.dtype or None, optional
            Floating dtype for internal computations.  Defaults to
            ``float32`` on CUDA and ``float64`` on CPU.

        Returns
        -------
        re, im : torch.Tensor
            Tensors with the same batch shape as the input parameters,
            containing the real and imaginary parts of the integral.

        USAGE
        =====
        1. 1D integral without parameters (single value)::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Integral, exp, oo

            x = symbols("x", real=True)
            integral_expr = Integral(exp(-x**2), (x, -oo, oo))

            lt = libtorch()
            texpr = lt.torchify(integral_expr)

            # No parameters → pass ``None`` for ``params_values``
            re, im = texpr.torch_integrate_batched(
                params_values=None,
                N=121,
            )
            ````

        2. 1D integral evaluated on a grid of parameter values::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Integral, exp, I, oo
            import torch

            x, k = symbols("x k", real=True)
            integral_expr = Integral(exp(-x**2) * exp(I * k * x), (x, -oo, oo))

            lt = libtorch()
            texpr = lt.torchify(integral_expr)

            # Parameter grid: 250 points between -5 and 5
            k_grid = torch.linspace(-5.0, 5.0, 250).unsqueeze(-1)  # shape (250, 1)

            re, im = texpr.torch_integrate_batched(
                params_values=k_grid,
                N=51,
                chunk_size_params=64,
                chunk_size_points=10_000,
            )
            ````
        """
        if self.dim <= 0:
            raise ValueError(f"self.dim must be positive; got dim={self.dim}")
        if self.n_params < 0:
            raise ValueError(f"self.n_params must be non-negative; got n_params={self.n_params}")
        if len(self.domain) != self.dim:
            raise ValueError(f"self.domain must have length dim={self.dim}; got {len(self.domain)}")

        if device is None:
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if dtype is None:
            dtype = torch.float32 if device.type == "cuda" else torch.float64

        if chunk_size_params <= 0:
            raise ValueError(f"chunk_size_params must be positive; got {chunk_size_params}")
        if chunk_size_points is not None and chunk_size_points <= 0:
            raise ValueError(f"chunk_size_points must be positive; got {chunk_size_points}")

        params_flat, batch_shape = self._normalize_params_values(
            params_values,
            self.n_params,
            device=device,
            dtype=dtype,
        )
        batch_size = int(params_flat.shape[0])

        rule = self._resolve_quadrature_rule(method)

        coords_1d = []
        weights_1d = []
        for (lower_bound, upper_bound) in self.domain:
            lower_bound = float(lower_bound)
            upper_bound = float(upper_bound)
            coords, weights = rule(lower_bound, upper_bound, N, device=device, dtype=dtype)
            coords = torch.as_tensor(coords, device=device, dtype=dtype).reshape(-1)
            weights = torch.as_tensor(weights, device=device, dtype=dtype).reshape(-1)
            if coords.numel() != weights.numel():
                raise ValueError(
                    f"Quadrature rule returned mismatched coordinates/weights lengths: {coords.numel()} vs {weights.numel()}"
                )
            coords_1d.append(coords)
            weights_1d.append(weights)

        if self.dim == 1:
            grid = coords_1d[0].unsqueeze(1)
            weights = weights_1d[0]
        else:
            grid = torch.cartesian_prod(*coords_1d)
            weights = weights_1d[0]
            for weight_1d in weights_1d[1:]:
                weights = torch.kron(weights, weight_1d)
        n_points = int(grid.shape[0])

        if chunk_size_points is None:
            chunk_size_points = n_points

        re_out = torch.empty(batch_size, device=device, dtype=dtype)
        im_out = torch.empty(batch_size, device=device, dtype=dtype)

        for start_param in range(0, batch_size, chunk_size_params):
            stop_param = min(batch_size, start_param + chunk_size_params)
            param_chunk = params_flat[start_param:stop_param, :]
            param_args = [param_chunk[:, index].unsqueeze(0) for index in range(self.n_params)]

            re_acc = torch.zeros(stop_param - start_param, device=device, dtype=dtype)
            im_acc = torch.zeros(stop_param - start_param, device=device, dtype=dtype)

            for start_point in range(0, n_points, chunk_size_points):
                stop_point = min(n_points, start_point + chunk_size_points)
                grid_chunk = grid[start_point:stop_point, :]
                weight_chunk = weights[start_point:stop_point].unsqueeze(1)

                var_args = [grid_chunk[:, index].unsqueeze(1) for index in range(self.dim)]
                values = self.func(*var_args, *param_args)
                values = torch.as_tensor(values, device=device)

                if values.ndim == 1:
                    values = values.unsqueeze(1)
                elif (
                    values.ndim == 2
                    and values.shape[0] == (stop_param - start_param)
                    and values.shape[1] == (stop_point - start_point)
                ):
                    values = values.t().contiguous()

                if values.shape[0] != (stop_point - start_point) or values.shape[1] != (stop_param - start_param):
                    raise RuntimeError(
                        f"Unexpected integrand output shape {tuple(values.shape)}; expected ({stop_point - start_point}, {stop_param - start_param})"
                    )

                if torch.is_complex(values):
                    re_acc += (weight_chunk * values.real).sum(dim=0)
                    im_acc += (weight_chunk * values.imag).sum(dim=0)
                else:
                    re_acc += (weight_chunk * values).sum(dim=0)

            re_out[start_param:stop_param] = re_acc
            im_out[start_param:stop_param] = im_acc

        re_out = re_out.reshape(batch_shape)
        im_out = im_out.reshape(batch_shape)
        return re_out, im_out

    def torch_integrate_batched_simpson(
        self,
        *,
        params_values = None,
        N: int = 121,
        chunk_size_params: int = 256,
        chunk_size_points: int | None = None,
        device=None,
        dtype=None,
    ):
        """Backward-compatible wrapper for ``torch_integrate_batched(method='simpson')``."""
        return self.torch_integrate_batched(
            params_values=params_values,
            method="simpson",
            N=N,
            chunk_size_params=chunk_size_params,
            chunk_size_points=chunk_size_points,
            device=device,
            dtype=dtype,
        )
    


class libtorch:
    def __init__(self):
        pass

    @staticmethod
    def _normalize_cov_method(change_of_variables_method: str) -> str:
        method = str(change_of_variables_method).strip().lower()
        aliases = {
            "tangent": "tangent",
            "tan": "tangent",
            "algebraic": "algebraic",
            "rational": "algebraic",
            "tanh-sinh": "tanh-sinh",
            "tanh_sinh": "tanh-sinh",
            "tanhsinh": "tanh-sinh",
        }
        if method not in aliases:
            raise ValueError(
                "Unknown change_of_variables_method. Expected one of {'tangent', 'algebraic', 'tanh-sinh'}. "
                f"Got: {change_of_variables_method}"
            )
        return aliases[method]

    @staticmethod
    def _inset_open_interval(lower: float, upper: float, eps: float) -> list[float]:
        return [float(lower + eps), float(upper - eps)]

    def _transform_limit_with_method(self, variable, lower, upper, *, change_of_variables_method: str, eps: float):
        """Return transformed variable, substitution, jacobian and finite domain."""
        method = self._normalize_cov_method(change_of_variables_method)

        if lower != -oo and upper != oo:
            return variable, variable, Integer(1), [float(lower), float(upper)]

        t_v = Symbol(f"t_{variable.name}", real=True)
        open_interval = None

        if method == "tangent":
            if lower == -oo and upper == oo:
                mapped = tan(t_v)
                open_interval = (float(-pi / 2), float(pi / 2))
            elif lower != -oo and upper == oo:
                mapped = lower + tan(t_v) ** 2
                open_interval = (0.0, float(pi / 2))
            elif lower == -oo and upper != oo:
                mapped = upper - tan(t_v) ** 2
                open_interval = (0.0, float(pi / 2))
            else:
                raise ValueError("Unexpected limit pattern in tangent transform")

        elif method == "algebraic":
            if lower == -oo and upper == oo:
                mapped = t_v / (1 - t_v ** 2)
                open_interval = (-1.0, 1.0)
            elif lower != -oo and upper == oo:
                mapped = lower + t_v / (1 - t_v)
                open_interval = (0.0, 1.0)
            elif lower == -oo and upper != oo:
                mapped = upper - t_v / (1 - t_v)
                open_interval = (0.0, 1.0)
            else:
                raise ValueError("Unexpected limit pattern in algebraic transform")

        elif method == "tanh-sinh":
            # Finite-interval parameter t in (-1,1) or (0,1) mapped through sinh(atanh(t)).
            core = sinh(atanh(t_v))
            if lower == -oo and upper == oo:
                mapped = core
                open_interval = (-1.0, 1.0)
            elif lower != -oo and upper == oo:
                mapped = lower + core ** 2
                open_interval = (0.0, 1.0)
            elif lower == -oo and upper != oo:
                mapped = upper - core ** 2
                open_interval = (0.0, 1.0)
            else:
                raise ValueError("Unexpected limit pattern in tanh-sinh transform")

        else:
            raise ValueError(f"Unsupported change_of_variables_method: {method}")

        jacobian = diff(mapped, t_v)
        domain = self._inset_open_interval(open_interval[0], open_interval[1], eps)
        return t_v, mapped, jacobian, domain

    def _merge_lambdify_modules(self, modules, extra_mapping):
        """Merge lambdify modules with an extra mapping dict at highest priority."""
        if not extra_mapping:
            return modules

        if modules is None:
            modules = self._default_modules()

        if isinstance(modules, (list, tuple)):
            merged = list(modules)
            merged.insert(0, dict(extra_mapping))
            return merged

        return [dict(extra_mapping), modules]

    def _build_nested_definite_integral_callable(self, texpr: TorchExpr, n_params: int, *, inner_N: int = 61):
        """Create a callable for lambdify that evaluates a definite nested integral numerically."""

        def eval_nested(*args):
            if len(args) != n_params:
                raise ValueError(f"Expected {n_params} nested integral params; got {len(args)}")

            if n_params == 0:
                re_val, im_val = texpr.torchquad_integrate(N=inner_N)
                re_tensor = torch.as_tensor(re_val)
                im_tensor = torch.as_tensor(im_val, device=re_tensor.device)
                return re_tensor + 1j * im_tensor

            args_tensors = [torch.as_tensor(arg) for arg in args]
            broadcasted = torch.broadcast_tensors(*args_tensors)
            target_device = broadcasted[0].device
            batch_shape = tuple(broadcasted[0].shape)
            flat_params = [tensor.reshape(-1) for tensor in broadcasted]
            batch_size = int(flat_params[0].shape[0])

            out_re = []
            out_im = []
            for idx in range(batch_size):
                params_i = [flat_params[param_index][idx] for param_index in range(n_params)]
                re_val, im_val = texpr.torchquad_integrate(params_values=params_i, N=inner_N)
                out_re.append(torch.as_tensor(re_val, device=target_device))
                out_im.append(torch.as_tensor(im_val, device=target_device))

            re_tensor = torch.stack(out_re).reshape(batch_shape)
            im_tensor = torch.stack(out_im).reshape(batch_shape).to(device=re_tensor.device)
            return re_tensor + 1j * im_tensor

        return eval_nested

    def _replace_nested_definite_integrals(
        self,
        expr,
        *,
        inner_N: int = 61,
        change_of_variables_method: str = "tangent",
        cov_eps: float = 1e-7,
    ):
        """Replace nested definite Integrals with numeric callables.

        Any nested integral with non-definite limits triggers a ValueError.
        """
        if not hasattr(expr, "has") or not expr.has(Integral):
            return expr, {}

        expr_work = expr
        nested_mapping = {}
        counter = 0

        while expr_work.has(Integral):
            leaf_integrals = [itg for itg in expr_work.atoms(Integral) if not itg.function.has(Integral)]
            if not leaf_integrals:
                break

            for inner_integral in leaf_integrals:
                for lim in inner_integral.limits:
                    if len(lim) != 3:
                        raise ValueError(
                            "The integrand must not contain unevaluated indefinite integrals. "
                            f"Found nested integral with non-definite limit: {inner_integral}"
                        )

                integration_vars = [lim[0] for lim in inner_integral.limits]
                param_symbols = sorted(
                    list(inner_integral.free_symbols - set(integration_vars)),
                    key=lambda symbol: symbol.name,
                )

                nested_texpr = self.torchify(
                    inner_integral,
                    change_of_variables_method=change_of_variables_method,
                    cov_eps=cov_eps,
                )
                nested_name = f"_nested_definite_integral_{counter}"
                counter += 1

                nested_mapping[nested_name] = self._build_nested_definite_integral_callable(
                    nested_texpr,
                    len(param_symbols),
                    inner_N=inner_N,
                )

                placeholder = Function(nested_name)(*param_symbols)
                expr_work = expr_work.xreplace({inner_integral: placeholder})

        if expr_work.has(Integral):
            raise ValueError(
                "The integrand must not contain unevaluated integrals. "
                "Only nested definite integrals are supported."
            )

        return expr_work, nested_mapping

    def _default_modules(self):
        """Internal helper: default SymPy → torch mapping for ``lambdify``.

        The returned list is passed directly to SymPy's ``lambdify`` as the
        ``modules`` argument.  It exposes a minimal subset of functions
        implemented with torch operations so that the resulting numerical
        function is differentiable and GPU-friendly where possible.

        Keeping this mapping in a separate method avoids cluttering
        :meth:`torchify` and makes it easy to customize or extend in user
        code by subclassing ``libtorch``.
        """

        def safe_sqrt(x):
            return torch.sqrt(torch.as_tensor(x, dtype=torch.float64))

        def safe_erf(x):
            """Error function that supports real torch tensors and complex inputs."""
            if torch.is_tensor(x):
                if torch.is_complex(x):
                    # torch.erf is not implemented for complex tensors on CPU.
                    out = scipy_special.erf(x.detach().cpu().numpy())
                    return torch.as_tensor(out, device=x.device, dtype=x.dtype)
                return torch.erf(x)
            return scipy_special.erf(x)

        def safe_erfc(x):
            """Complementary error function with complex-input fallback."""
            if torch.is_tensor(x):
                if torch.is_complex(x):
                    out = scipy_special.erfc(x.detach().cpu().numpy())
                    return torch.as_tensor(out, device=x.device, dtype=x.dtype)
                return torch.erfc(x)
            return scipy_special.erfc(x)

        return [{
            # Trigonometric
            "sin": torch.sin, "cos": torch.cos, "tan": torch.tan,
            "asin": torch.asin, "acos": torch.acos, "atan": torch.atan,
            "atan2": torch.atan2,
            # Hyperbolic
            "sinh": torch.sinh, "cosh": torch.cosh, "tanh": torch.tanh,
            "asinh": torch.asinh, "acosh": torch.acosh, "atanh": torch.atanh,
            # Exponentials / logs
            "exp": torch.exp, "log": torch.log, "ln": torch.log,
            "log10": torch.log10, "log2": torch.log2,
            # Special functions used by symbolic integral evaluation
            "erf": safe_erf, "erfc": safe_erfc,
            # Roots / powers
            "sqrt": safe_sqrt, "Pow": torch.pow,
            # Misc
            "Abs": torch.abs, "sign": torch.sign,
            "conjugate": torch.conj, "conj": torch.conj,
            "floor": torch.floor, "ceiling": torch.ceil,
            "Min": torch.minimum, "Max": torch.maximum,
            # Piecewise / heaviside
            "Heaviside": lambda x: torch.heaviside(x, torch.zeros_like(x)),
            # Constants
            "pi": getattr(torch, "pi", float(np.pi)),
            "E": float(np.e),
        }]

    def torchify(
        self,
        expr,
        *,
        variables=None,
        limits=None,
        params=None,
        modules=None,
        change_of_variables_method: str = "tangent",
        cov_eps: float = 1e-7,
    ):
        """Convert a SymPy object into a torch-compatible callable or TorchExpr.

        This mirrors the original functional ``torchify`` but lives as a
        method on ``libtorch``.  It performs change-of-variables for
        infinite / semi-infinite limits at the SymPy level (once), then
        ``lambdify``'s the transformed expression.

        ``expr`` can be:

        - a plain SymPy expression (e.g. ``exp(-x**2)``),
        - a SymPy ``Integral`` (``Integral(f(x), (x, -oo, oo))``), or
        - a SymPy ``Eq`` whose right-hand side is an ``Integral``.

        If an ``Integral`` (or equation with an integral) is passed, the
        integrand and limits are extracted internally, so you no longer
        need to manually use ``.function`` / ``.limits``.

        Parameters
        ----------
        expr : sympy.Expr | sympy.Integral | sympy.Eq
            Symbolic expression (may be complex-valued) or an integral.
        variables : list[sympy.Symbol] or None
            Integration variables (order matters).  If ``expr`` is an
            ``Integral`` and ``variables`` is ``None``, they are inferred
            from the integral limits.
        limits : list[tuple] or None
            ``(var, lower, upper)`` for each variable.  ``oo`` / ``-oo``
            trigger automatic change of variables.  If ``expr`` is an
            ``Integral`` and ``limits`` is ``None``, they are taken from
            ``expr.limits``.
        params : list[sympy.Symbol] or None
            Extra symbolic parameters that appear in the transformed
            expression but are **not** integrated over.  If ``None``, they
            are inferred automatically from the free symbols.
        modules : list[dict] or None
            Custom ``lambdify`` module list.  ``None`` → default torch
            mapping from ``_default_modules``.
        change_of_variables_method : str, optional
            Method used to transform infinite/semi-infinite domains to finite
            ones. Supported values: ``"tangent"``, ``"algebraic"``,
            ``"tanh-sinh"``.
        cov_eps : float, optional
            Small inward shift for transformed open intervals to avoid
            endpoint singular evaluations.

        Returns
        -------
        TorchExpr
            If limits are provided (directly or via an ``Integral``), ready
            for numerical integration.
        callable
            If ``limits`` is ``None``, a plain torch-ready function
            (no integration).

        USAGE
        =====
        1. Plain function (no integration)::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, exp

            x = symbols("x", real=True)

            lt = libtorch()
            f  = lt.torchify(exp(-x**2), variables=[x])
            ````

        2. Integral with finite limits::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Integral, exp

            x = symbols("x", real=True)
            integral = Integral(exp(-x**2), (x, 0, 1))

            lt    = libtorch()
            texpr = lt.torchify(integral)
            ````

        3. Integral with infinite limits (automatic change of variables)::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Integral, exp, oo

            x = symbols("x", real=True)
            integral = Integral(exp(-x**2), (x, -oo, oo))

            lt    = libtorch()
            texpr = lt.torchify(integral)
            ````

        4. Equation with integral on the right-hand side::

            ````python
            from libphysics.libtorch import libtorch
            from sympy import symbols, Eq, Integral, exp, oo

            x, y = symbols("x y", real=True)
            eq = Eq(y, Integral(exp(-x**2), (x, -oo, oo)))

            lt    = libtorch()
            texpr = lt.torchify(eq)  # integral auto-extracted
            ````
        """

        # Detect integral / equation-with-integral and extract base expr/limits.
        integral = None
        if isinstance(expr, Integral):
            integral = expr
        elif isinstance(expr, Eq) and isinstance(expr.rhs, Integral):
            integral = expr.rhs

        if integral is not None:
            base_expr = integral.function
            nested_modules = {}
            if base_expr.has(Integral):
                # Evaluate nested definite integrals numerically using torchquad.
                base_expr, nested_modules = self._replace_nested_definite_integrals(
                    base_expr,
                    change_of_variables_method=change_of_variables_method,
                    cov_eps=cov_eps,
                )
            if limits is None:
                limits = list(integral.limits)
            if variables is None and limits is not None:
                variables = [lim[0] for lim in limits]
        else:
            base_expr = expr
            nested_modules = {}

        if variables is None:
            raise ValueError("'variables' must be provided when expr is not an Integral")

        if modules is None:
            modules = self._default_modules()
        modules = self._merge_lambdify_modules(modules, nested_modules)

        variables = list(variables)

        # ------------------------------------------------------------------
        # Simple case: no limits → plain lambdify
        # ------------------------------------------------------------------
        if limits is None:
            return lambdify(tuple(variables), base_expr, modules=modules)

        # ------------------------------------------------------------------
        # With limits: symbolic change of variables (done once, fast forever)
        # ------------------------------------------------------------------
        limit_map = {lim[0]: (lim[1], lim[2]) for lim in limits}

        new_vars = []
        domain = []
        subs_map = {}
        jacobian = Integer(1)

        for v in variables:
            lower, upper = limit_map[v]
            t_var, mapped_expr, jac_expr, domain_entry = self._transform_limit_with_method(
                v,
                lower,
                upper,
                change_of_variables_method=change_of_variables_method,
                eps=cov_eps,
            )
            new_vars.append(t_var)
            domain.append(domain_entry)
            if t_var != v:
                subs_map[v] = mapped_expr
            jacobian *= jac_expr

        expr_work = base_expr.subs(subs_map) * jacobian

        if params is None:
            params = sorted(
                list(expr_work.free_symbols - set(new_vars)),
                key=lambda s: s.name,
            )
        else:
            params = list(params)

        # Lambdify: new_vars first, then params
        arglist = tuple(new_vars + params)
        func = lambdify(arglist, expr_work, modules=modules)

        return TorchExpr(
            func=func,
            domain=domain,
            dim=len(new_vars),
            n_params=len(params),
            sympy_expr=expr_work,
            variables=new_vars,
        )


# ---------------------------------------------------------------------------
# Functional compatibility helpers
# ---------------------------------------------------------------------------
def torchify(
    expr,
    variables=None,
    limits=None,
    params=None,
    modules=None,
    change_of_variables_method: str = "tangent",
    cov_eps: float = 1e-7,
):
    """Convert a SymPy object into a torch-compatible callable or ``TorchExpr``.

    This functional wrapper is kept for backward compatibility with the
    original ``libtorch`` API.  It delegates to ``libtorch().torchify(...)``
    so the functional and object-oriented styles share the same core logic.

    Parameters
    ----------
    expr : sympy.Expr | sympy.Integral | sympy.Eq
        Symbolic expression or integral-like SymPy object.
    variables : list[sympy.Symbol] or None, optional
        Variables passed through to ``libtorch.torchify``.
    limits : list[tuple] or None, optional
        Integration limits passed through to ``libtorch.torchify``.
    params : list[sympy.Symbol] or None, optional
        Additional symbolic parameters.
    modules : list[dict] or None, optional
        Custom ``lambdify`` modules mapping.

    Returns
    -------
    TorchExpr | callable
        Same return contract as ``libtorch().torchify(...)``.

    USAGE
    =====
    1. Functional style kept in parallel with the OOP API::

        ````python
        from libphysics.libtorch import torchify
        from sympy import symbols, Integral, exp, oo

        x = symbols("x", real=True)
        integral_expr = Integral(exp(-x**2), (x, -oo, oo))

        texpr = torchify(integral_expr)
        ````
    """
    return libtorch().torchify(
        expr,
        variables=variables,
        limits=limits,
        params=params,
        modules=modules,
        change_of_variables_method=change_of_variables_method,
        cov_eps=cov_eps,
    )


def torchquad_integrate(texpr: TorchExpr, params_values=None, method=None, N: int = 21):
    """Numerically integrate a ``TorchExpr`` using torchquad.

    This functional wrapper is kept for backward compatibility with the
    original ``libtorch`` API.  New code can call
    :meth:`TorchExpr.torchquad_integrate` directly.

    Parameters
    ----------
    texpr : TorchExpr
        Object returned by ``libtorch().torchify(..., limits=...)``.
    params_values : list | tuple | torch.Tensor | None, optional
        Numerical parameter values in the same order expected by
        ``texpr.func``.
    method : torchquad integrator instance or None, optional
        If ``None``, a default ``Simpson`` integrator is created.
    N : int, optional
        Resolution parameter passed to torchquad.

    Returns
    -------
    re, im : torch.Tensor
        Real and imaginary parts of the integral.

    USAGE
    =====
    1. Functional style kept in parallel with the OOP API::

        ````python
        from libphysics.libtorch import libtorch, torchquad_integrate
        from sympy import symbols, Integral, exp, oo

        x = symbols("x", real=True)
        integral_expr = Integral(exp(-x**2), (x, -oo, oo))

        lt = libtorch()
        texpr = lt.torchify(integral_expr)
        re, im = torchquad_integrate(texpr, N=121)
        ````
    """
    return texpr.torchquad_integrate(params_values=params_values, method=method, N=N)


def _simpson_weights_1d(a: float, b: float, N: int, *, device=None, dtype=None) -> torch.Tensor:
    """Backward-compatible Simpson weights wrapper."""
    _coords, weights = TorchExpr._rule_simpson(a, b, N, device=device, dtype=dtype)
    return weights


def _is_scalar_like(x) -> bool:
    """Functional wrapper around ``TorchExpr._is_scalar_like``."""
    return TorchExpr._is_scalar_like(x)


def _normalize_params_values(params_values, n_params: int, *, device, dtype):
    """Functional wrapper around ``TorchExpr._normalize_params_values``."""
    return TorchExpr._normalize_params_values(params_values, n_params, device=device, dtype=dtype)


def torch_integrate_batched(
    texpr: TorchExpr,
    params_values=None,
    *,
    method: str | Callable = "simpson",
    N: int = 121,
    chunk_size_params: int = 256,
    chunk_size_points: int | None = None,
    device=None,
    dtype=None,
):
    """Batched tensor-product quadrature integration for any dimension.

    This functional wrapper is kept so existing code that used the
    original function-based API continues to work.  New code can call
    :meth:`TorchExpr.torch_integrate_batched` directly.

    Parameters
    ----------
    texpr : TorchExpr
        Object returned by ``libtorch().torchify(..., limits=...)``.
    params_values : tensor | list | tuple | None
        Numerical parameter values.
    N : int, optional
        Odd number of Simpson points per dimension.
    chunk_size_params : int, optional
        Number of parameter points processed per chunk.
    chunk_size_points : int or None, optional
        Number of sample points processed per chunk.
    device : torch.device or None, optional
        Device used for internal tensors.
    dtype : torch.dtype or None, optional
        Floating dtype used for internal tensors.

    Returns
    -------
    re, im : torch.Tensor
        Real and imaginary parts of the integral.

    USAGE
    =====
    1. Functional style kept in parallel with the OOP API::

        ````python
        from libphysics.libtorch import libtorch, torch_integrate_batched_simpson
        from sympy import symbols, Integral, exp, oo

        x = symbols("x", real=True)
        integral_expr = Integral(exp(-x**2), (x, -oo, oo))

        lt = libtorch()
        texpr = lt.torchify(integral_expr)
        re, im = torch_integrate_batched_simpson(texpr, params_values=None, N=121)
        ````
    """
    return texpr.torch_integrate_batched(
        params_values=params_values,
        method=method,
        N=N,
        chunk_size_params=chunk_size_params,
        chunk_size_points=chunk_size_points,
        device=device,
        dtype=dtype,
    )


def torch_integrate_batched_simpson(
    texpr: TorchExpr,
    params_values=None,
    *,
    N: int = 121,
    chunk_size_params: int = 256,
    chunk_size_points: int | None = None,
    device=None,
    dtype=None,
):
    """Backward-compatible wrapper for ``torch_integrate_batched(method='simpson')``."""
    return texpr.torch_integrate_batched(
        params_values=params_values,
        method="simpson",
        N=N,
        chunk_size_params=chunk_size_params,
        chunk_size_points=chunk_size_points,
        device=device,
        dtype=dtype,
    )

lt = libtorch()