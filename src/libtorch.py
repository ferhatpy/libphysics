# -*- coding: utf-8 -*-
"""
libtorch.py — SymPy-to-Torch conversion and torchquad integration.

Part of libphysics.

Two-function design for performance:
    1. torchify()            — expensive SymPy work done ONCE (change-of-vars + lambdify)
    2. torchquad_integrate() — cheap torch-only integration, called MANY times

Usage
=====
        from libphysics.libtorch import torchify, torchquad_integrate

        # Step 1: build once (slow: SymPy subs + lambdify, ~200-500 ms)
        wigner = torchify(
                integrand,
                variables=[y_A, y_B],
                limits=[(y_A, -oo, oo), (y_B, -oo, oo)],
                params=[x_A, p_A, x_B, p_B],
        )

        # Step 2: integrate many times (fast: pure torch, ~3-8 ms each)
    re, im = torchquad_integrate(wigner, params_values=[1, 1, 1, 1], N=121)
"""
from dataclasses import dataclass
import sys
from typing import List, Callable

import numpy as np
import torch
from sympy import lambdify, Symbol, tan, cos, pi, oo, Integer

from loguru import logger

# Remove loguru's default handler only if user requests it via env var.
# Set LIBPHYSICS_REMOVE_LOGURU=1 or "true"/"yes" to disable.
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

# ---------------------------------------------------------------------------
# TorchExpr — container returned by torchify when limits are provided
# ---------------------------------------------------------------------------
@dataclass
class TorchExpr:
    """Pre-built torch function + finite integration domain."""
    func: Callable                # f(new_var_0, …, new_var_n, param_0, …, param_m)
    domain: List[List[float]]     # finite box for torchquad
    dim: int                      # number of integration variables
    n_params: int                 # number of extra parameters


# ---------------------------------------------------------------------------
# torchify  —  expensive work done ONCE
# ---------------------------------------------------------------------------
def torchify(expr, variables, limits=None, params=None, modules=None):
    """
    Convert a SymPy expression into a torch-compatible function via lambdify.

    If *limits* are provided, performs change-of-variables for infinite /
    semi-infinite limits at the **SymPy level** (symbolic, done once) and
    multiplies by the analytic Jacobian.  Returns a ``TorchExpr`` ready
    for ``torchquad_integrate``.

    If *limits* are **not** provided, returns a plain callable (simple
    lambdify with torch modules).

    Parameters
    ----------
    expr : sympy.Expr
        Symbolic expression (may be complex-valued).
    variables : list[sympy.Symbol]
        Integration variables (order matters).
    limits : list[tuple] or None
        ``(var, lower, upper)`` for each variable.  ``sympy.oo`` /
        ``-sympy.oo`` trigger automatic change of variables.
    params : list[sympy.Symbol] or None
        Extra symbolic parameters that appear in *expr* but are **not**
        integrated over.  If ``None``, inferred automatically.
    modules : list[dict] or None
        Custom lambdify module list.  ``None`` → default torch mapping.

    Returns
    -------
    TorchExpr   if *limits* were given  (use with ``torchquad_integrate``)
    callable    if *limits* were ``None`` (plain torch function)
    -------
    Usage: 
        x, k = sp.symbols("x k", real=True)
        integrand = sp.exp(-x**2) * sp.exp(sp.I * k * x)

        texpr = torchify(
            integrand,
            variables=[x],
            limits=[(x, -sp.oo, sp.oo)],
            params=[k],
            )

        # if the expr is given as an sp integral:
        inntegrand = integral_expression.function
        limits = list(integral_expression.limits)

        texpr = torchify(
            inntegrand,
            variables=[x],
            limits=limits,
            params=[k],
            )

        # if the expr is given as an sp equation integral:
        inntegrand = integral_equation.rhs.function
        limits = list(integral_equation.rhs.limits)

        texpr = torchify(
            inntegrand,
            variables=[x],
            limits=limits,
            params=[k],
            )
    """
    def safe_sqrt(x):
        return torch.sqrt(torch.as_tensor(x, dtype=torch.float64))
    if modules is None:
        modules = [{
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
            # Roots / powers
            "sqrt": safe_sqrt, "Pow": torch.pow,
            # Misc
            "Abs": torch.abs, "sign": torch.sign,
            "floor": torch.floor, "ceiling": torch.ceil,
            "Min": torch.minimum, "Max": torch.maximum,
            # Piecewise / heaviside
            "Heaviside": lambda x: torch.heaviside(x, torch.zeros_like(x)),
            # Constants
            "pi": getattr(torch, "pi", float(np.pi)),
            "E": float(np.e),
        }]

    variables = list(variables)

    # ------------------------------------------------------------------
    # Simple case: no limits → plain lambdify
    # ------------------------------------------------------------------
    if limits is None:
        return lambdify(tuple(variables), expr, modules=modules)

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

        if lower == -oo and upper == oo:
            t_v = Symbol(f"t_{v.name}", real=True)
            subs_map[v] = tan(t_v)
            jacobian *= 1 / cos(t_v) ** 2
            new_vars.append(t_v)
            domain.append([float(-pi / 2), float(pi / 2)])

        elif lower != -oo and upper == oo:
            t_v = Symbol(f"t_{v.name}", real=True)
            subs_map[v] = lower + tan(t_v) ** 2
            jacobian *= 2 * tan(t_v) / cos(t_v) ** 2
            new_vars.append(t_v)
            domain.append([0.0, float(pi / 2)])

        elif lower == -oo and upper != oo:
            t_v = Symbol(f"t_{v.name}", real=True)
            subs_map[v] = upper - tan(t_v) ** 2
            jacobian *= -2 * tan(t_v) / cos(t_v) ** 2
            new_vars.append(t_v)
            domain.append([0.0, float(pi / 2)])

        else:
            # finite [a, b]
            new_vars.append(v)
            domain.append([float(lower), float(upper)])

    expr_work = expr.subs(subs_map) * jacobian

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

    return TorchExpr(func=func, domain=domain, dim=len(new_vars), n_params=len(params))


# ---------------------------------------------------------------------------
# torchquad_integrate  —  cheap work, called MANY times
# ---------------------------------------------------------------------------
def torchquad_integrate(texpr, params_values=None, method=None, N=21):
    """
    Numerically integrate a ``TorchExpr`` using torchquad.

    Parameters
    ----------
    texpr : TorchExpr
        Object returned by ``torchify(..., limits=...)``.
    params_values : list or None
        Numerical values for the parameters (same order as *params* in torchify).
    method : torchquad integrator or None
        Defaults to ``Simpson()``.
    N : int
        Integrator resolution (torchquad's *N* parameter).

    Returns
    -------
    re, im : torch.Tensor
        Real and imaginary parts of the integral.
    -------
    Usage:
        re, im = torchquad_integrate(texpr, params_values=[1, 1, 1, 1], N=121)
    """
    from torchquad import Simpson

    if method is None:
        method = Simpson()

    param_vals = list(params_values) if params_values else []


    def f(d):
        vals = texpr.func(*[d[:, i] for i in range(d.shape[1])], *param_vals)
        return torch.as_tensor(vals, device=d.device)
            
    res = method.integrate(f, dim=texpr.dim, N=N, integration_domain= texpr.domain)
    return (res.real if torch.is_complex(res) else res), (res.imag if torch.is_complex(res) else torch.zeros_like(res))



def _simpson_weights_1d(a: float, b: float, N: int, *, device=None, dtype=None) -> torch.Tensor:
    """Return Simpson weights including the dx/3 scaling (shape: [N]).
     this is used to build the full tensor-product weights for torchquad_integrate.
    """
    if N < 3 or (N % 2) == 0:
        raise ValueError(f"Simpson rule requires odd N>=3; got N={N}")
    dx = (b - a) / (N - 1)
    w = torch.ones(N, device=device, dtype=dtype)
    # 1, 4, 2, 4, ..., 2, 4, 1 this is the simppson pattern
    w[1:-1:2] = 4
    w[2:-1:2] = 2
    w = w * (dx / 3.0)
    return w


def _is_scalar_like(x) -> bool:
    """True for Python numbers / 0-d tensors."""
    if torch.is_tensor(x):
        return x.ndim == 0
    return isinstance(x, (int, float, complex, np.number))


def _normalize_params_values(params_values, n_params: int, *, device, dtype):
    """Normalize params into (B, n_params) plus batch_shape.

    Accepts:
      - tensor of shape (..., n_params)
      - list/tuple length n_params containing tensors and/or scalars

    Scalars are broadcast to the inferred batch_shape.
    """
    if params_values is None:
        if n_params != 0:
            raise ValueError(f"Expected {n_params} params; got None")
        return torch.empty((1, 0), device=device, dtype=dtype), tuple()

    # Case A: tensor (..., n_params)
    if torch.is_tensor(params_values):
        params_tensor = params_values
        if params_tensor.shape[-1] != n_params:
            raise ValueError(
                f"params_values last dimension must be n_params={n_params}; got shape={tuple(params_tensor.shape)}"
            )
        batch_shape = tuple(params_tensor.shape[:-1])
        B = int(np.prod(batch_shape)) if batch_shape else 1
        params_flat = params_tensor.reshape(B, n_params).to(device=device, dtype=dtype)
        return params_flat, batch_shape

    # Case B: list/tuple of params (scalars and/or tensors)
    if not isinstance(params_values, (list, tuple)):
        # single scalar treated as 1 param
        params_values = [params_values]

    if len(params_values) != n_params:
        raise ValueError(f"Expected {n_params} params; got {len(params_values)}")

    # Infer batch_shape from the first non-scalar param (or scalars-only => ())
    batch_shape = None
    for p in params_values:
        if torch.is_tensor(p) and p.ndim > 0:
            batch_shape = tuple(p.shape)
            break
        if not _is_scalar_like(p):
            pt = torch.as_tensor(p)
            if pt.ndim > 0:
                batch_shape = tuple(pt.shape)
                break
    if batch_shape is None:
        batch_shape = tuple()

    # Broadcast all params to batch_shape, then stack
    params_list = []
    for p in params_values:
        pt = torch.as_tensor(p, device=device)
        if pt.ndim == 0:
            if batch_shape:
                pt = pt.expand(batch_shape)
        else:
            if tuple(pt.shape) != batch_shape:
                raise ValueError(
                    f"All parameter tensors must have the same shape; got {tuple(pt.shape)} vs {batch_shape}"
                )
        params_list.append(pt)

    params_tensor = torch.stack(params_list, dim=-1).to(dtype=dtype)
    B = int(np.prod(batch_shape)) if batch_shape else 1
    params_flat = params_tensor.reshape(B, n_params)
    return params_flat, batch_shape


def torch_integrate_batched_simpson(
    texpr: TorchExpr,
    params_values,
    *,
    N: int = 121,
    chunk_size_params: int = 256,
    chunk_size_points: int | None = None,
    device=None,
    dtype=None,
):
    """Batched tensor-product Simpson integration for any dimension.

    This integrates a ``TorchExpr`` over its finite domain using a tensor-product
    Simpson rule with *N* points per dimension.

    Supports a batch of parameter points (e.g. a 250x250 mesh) by keeping the
    batch dimension and summing only over the sample grid.

    Notes
    -----
    - The number of sample points grows as $N^{dim}$. For dim>3 this can become
      very large quickly. Use smaller N and/or set ``chunk_size_points``.
    - ``chunk_size_points`` trades more Python overhead for lower peak memory.

    Parameters
    ----------
    texpr : TorchExpr
        Output of ``torchify(..., limits=...)``.
    params_values : tensor | list/tuple
        Either a tensor of shape ``(..., n_params)`` or a list/tuple of length
        ``n_params`` of tensors/scalars with broadcastable batch shapes.
    N : int
        Odd Simpson points per dimension.
    chunk_size_params : int
        Number of parameter points processed per chunk.
    chunk_size_points : int | None
        Optional number of sample points processed per chunk.

    Returns
    -------
    re, im : torch.Tensor
        Shapes match the parameter batch shape (e.g. 250x250).
    -------
    Usage:
        >>> re, im = torch_integrate_batched_simpson(
        ...     texpr,
        ...     params_values=params_grid,  # shape (250, 250, n_params)
        ...     N=21,
        ...     chunk_size_params=256,
        ...      chunk_size_points=10000)
    """
    if texpr.dim <= 0:
        raise ValueError(f"texpr.dim must be positive; got dim={texpr.dim}")
    if texpr.n_params < 0:
        raise ValueError(f"texpr.n_params must be non-negative; got n_params={texpr.n_params}")
    if len(texpr.domain) != texpr.dim:
        raise ValueError(f"texpr.domain must have length dim={texpr.dim}; got {len(texpr.domain)}")

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if dtype is None:
        dtype = torch.float32 if device.type == "cuda" else torch.float64

    if chunk_size_params <= 0:
        raise ValueError(f"chunk_size_params must be positive; got {chunk_size_params}")
    if chunk_size_points is not None and chunk_size_points <= 0:
        raise ValueError(f"chunk_size_points must be positive; got {chunk_size_points}")

    params_flat, batch_shape = _normalize_params_values(
        params_values,
        texpr.n_params,
        device=device,
        dtype=dtype,
    )
    B = int(params_flat.shape[0])

    # Build per-dimension coordinates and weights
    coords_1d = []
    weights_1d = []
    for (a, b) in texpr.domain:
        a = float(a)
        b = float(b)
        coords_1d.append(torch.linspace(a, b, N, device=device, dtype=dtype))
        weights_1d.append(_simpson_weights_1d(a, b, N, device=device, dtype=dtype))

    # Build full grid points (P, dim) and weights (P,)
    # This is the most straightforward approach; chunk_size_points can reduce peak memory.
    grid = torch.cartesian_prod(*coords_1d)  # (P, dim)
    w = weights_1d[0]
    for wi in weights_1d[1:]:
        w = torch.kron(w, wi)
    P = int(grid.shape[0])

    if chunk_size_points is None:
        chunk_size_points = P

    re_out = torch.empty(B, device=device, dtype=dtype)
    im_out = torch.empty(B, device=device, dtype=dtype)

    # Evaluate in chunks over params and points
    for start_p in range(0, B, chunk_size_params):
        stop_p = min(B, start_p + chunk_size_params)
        pc = params_flat[start_p:stop_p, :]  # (Bc, n_params)
        param_args = [pc[:, i].unsqueeze(0) for i in range(texpr.n_params)]

        re_acc = torch.zeros(stop_p - start_p, device=device, dtype=dtype)
        im_acc = torch.zeros(stop_p - start_p, device=device, dtype=dtype)

        for start_x in range(0, P, chunk_size_points):
            stop_x = min(P, start_x + chunk_size_points)
            g = grid[start_x:stop_x, :]  # (Px, dim)
            ww = w[start_x:stop_x].unsqueeze(1)  # (Px, 1)

            var_args = [g[:, i].unsqueeze(1) for i in range(texpr.dim)]
            vals = texpr.func(*var_args, *param_args)
            vals = torch.as_tensor(vals, device=device)

            # Expect (Px, Bc) after broadcasting.
            if vals.ndim == 1:
                vals = vals.unsqueeze(1)
            elif vals.ndim == 2 and vals.shape[0] == (stop_p - start_p) and vals.shape[1] == (stop_x - start_x):
                # Sometimes lambdify might return (Bc, Px)
                vals = vals.t().contiguous()

            if vals.shape[0] != (stop_x - start_x) or vals.shape[1] != (stop_p - start_p):
                raise RuntimeError(
                    f"Unexpected integrand output shape {tuple(vals.shape)}; expected ({stop_x - start_x}, {stop_p - start_p})"
                )

            if torch.is_complex(vals):
                re_acc += (ww * vals.real).sum(dim=0)
                im_acc += (ww * vals.imag).sum(dim=0)
            else:
                re_acc += (ww * vals).sum(dim=0)
                # im stays zero

        re_out[start_p:stop_p] = re_acc
        im_out[start_p:stop_p] = im_acc

    re_out = re_out.reshape(batch_shape)
    im_out = im_out.reshape(batch_shape)
    return re_out, im_out

