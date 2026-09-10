# -*- coding: utf-8 -*-
"""torchsympy_new.py — SymPy-to-Torch conversion and numerical integration helpers.

Standalone library (no package, a single importable module).

WHAT IS NEW IN THIS FILE
========================
This is torchsympy.py with one addition: :meth:`TorchSymPy.eval_numeric` now
recognises *separable* integrands and evaluates them as a product of
one-dimensional integrals, each refined until it converges, instead of using a
single tensor-product rule with a fixed node budget shared between the axes.

That is what makes oscillatory diffraction integrands come out right.  A 2D
Fresnel integral evaluated with ``N = 5001`` really gets only ~70 nodes per
axis, far fewer than the number of fringes of the integrand, so the quadrature
aliased them into a spurious pattern along the four sides of the screen at
z = 400 mm, and into noise at z = 40 mm.  The rewritten path uses the full node
count on each axis and picks it automatically, giving ~1e-13 relative accuracy
at every distance (see ``fresnel_torchsympy_check.py``).

The public API, the solvers and every other code path are unchanged;
non-separable integrands (e.g. the Rayleigh-Sommerfeld kernel) still go through
the original solvers, and ``TorchSymPy(separable_integration=False)`` restores
the previous behaviour exactly.  See the long comment above
:func:`_separable_plan` for the details.

Two-stage design for performance:
    1. ``TorchSymPy().torchify()``          — expensive SymPy work done ONCE
                                              (change-of-variables + ``lambdify``)
    2. ``TorchExpr.torch_integrate_batched()`` /
       ``TorchExpr.torchquad_integrate()``  — cheap torch-only integration,
                                              called MANY times

USAGE
=====
    import torchsympy

    # Step 1: build once (slow: SymPy subs + lambdify)
    lt = torchsympy.TorchSymPy()
    texpr = lt.torchify(
        sp.Integral(integrand, (y_A, -sp.oo, sp.oo), (y_B, -sp.oo, sp.oo))
    )

    # Step 2: integrate many times (fast: pure torch)
    re, im = texpr.torchquad_integrate(N=121)

    # Batched integration over a parameter grid (tensor-product quadrature)
    re_b, im_b = texpr.torch_integrate_batched(
        params_values=params_grid,
        N=21,
        chunk_size_params=256,
        chunk_size_points=10_000,
    )

Caching
=======
Quadrature nodes/weights, tensor-product grids and compiled (``lambdify``'d)
expressions are memoised, so repeated integrations with the same rule, domain
and resolution skip all the set-up work.  Call :func:`clear_caches` to release
that memory (useful on the GPU, or between benchmark runs).

See ``IMPROVEMENTS.md`` for the list of bugs fixed and speed-ups made in this
revision, and ``tests/test_torchsympy.py`` for the corresponding regression
tests.
"""
from __future__ import annotations

import math
from collections import OrderedDict
from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations, product
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import torch
from loguru import logger
from scipy import special as scipy_special
from sympy import (
    Add,
    Dummy,
    Eq,
    Float as SympyFloat,
    Function,
    I as sympy_I,
    Integer,
    Integral,
    Mul,
    Symbol,
    atanh,
    conjugate as sympy_conjugate,
    cos as sympy_cos,
    count_ops,
    diff,
    exp as sympy_exp,
    expand as sympy_expand,
    im as sympy_im,
    lambdify,
    oo,
    pi,
    re as sympy_re,
    simplify as sympy_simplify,
    sin as sympy_sin,
    sinh,
    tan,
)

__all__ = [
    "MAX_ELEMENTS_PER_EVALUATION",
    "TorchExpr",
    "TorchSymPy",
    "clear_caches",
    "setup_logging",
]

QuadratureRule = Callable[..., Tuple[torch.Tensor, torch.Tensor]]


# ----------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------
_LOG_HANDLER_ID: Optional[int] = None


def setup_logging(enable: bool = False, level: str = "INFO", *, sink=None) -> None:
    """Enable or disable the log messages emitted by this module.

    The library is silent by default.  Unlike a global ``logger.remove()`` this
    only touches this module's own records, so importing ``torchsympy`` never
    reconfigures the logging of the host application.

    Parameters
    ----------
    enable : bool
        ``True`` activates this module's log records.
    level : str
        Minimum level for the handler installed by this function.
    sink : file-like or None
        ``None`` (the default) re-enables the records and lets the loguru
        handlers configured by the application deal with them.  Passing a sink
        (e.g. ``sys.stderr``) additionally installs one dedicated handler,
        filtered to this module and limited to ``level``; it is removed again
        by the next call.

    Examples
    --------
    >>> setup_logging()                                        # silent (default)
    >>> setup_logging(enable=True)                             # use the app's handlers
    >>> setup_logging(enable=True, level="DEBUG", sink=sys.stderr)
    """
    global _LOG_HANDLER_ID

    if _LOG_HANDLER_ID is not None:
        try:
            logger.remove(_LOG_HANDLER_ID)
        except ValueError:  # pragma: no cover - handler already gone
            pass
        _LOG_HANDLER_ID = None

    if not enable:
        logger.disable(__name__)
        return

    logger.enable(__name__)
    if sink is not None:
        _LOG_HANDLER_ID = logger.add(
            sink, level=level, filter=lambda record: record["name"] == __name__
        )


setup_logging(enable=False)


# ----------------------------------------------------------------------
# dtype helpers
# ----------------------------------------------------------------------
_COMPLEX_OF_REAL: Dict[torch.dtype, torch.dtype] = {
    torch.float16: torch.complex32,
    torch.bfloat16: torch.complex64,
    torch.float32: torch.complex64,
    torch.float64: torch.complex128,
}
_REAL_OF_COMPLEX: Dict[torch.dtype, torch.dtype] = {
    torch.complex32: torch.float16,
    torch.complex64: torch.float32,
    torch.complex128: torch.float64,
}


def _real_dtype(dtype: Optional[torch.dtype]) -> Optional[torch.dtype]:
    """Real counterpart of ``dtype`` (identity for real dtypes)."""
    return _REAL_OF_COMPLEX.get(dtype, dtype)


def _complex_dtype(dtype: Optional[torch.dtype]) -> Optional[torch.dtype]:
    """Complex counterpart of ``dtype`` (identity for complex dtypes)."""
    if dtype is None:
        return None
    return _COMPLEX_OF_REAL.get(dtype, dtype)


def _as_float_tensor(x: Any) -> torch.Tensor:
    """Tensor view of ``x`` promoted to a floating/complex dtype if needed."""
    tensor = x if torch.is_tensor(x) else torch.as_tensor(x)
    if tensor.is_floating_point() or torch.is_complex(tensor):
        return tensor
    return tensor.to(torch.float64)


# ----------------------------------------------------------------------
# torch implementations of the SymPy functions used by ``lambdify``
# ----------------------------------------------------------------------
def _torch_sqrt(x: Any) -> torch.Tensor:
    """Square root that preserves dtype, device and complex values.

    The previous implementation forced ``dtype=torch.float64``, which discarded
    the imaginary part of complex integrands and upcast float32/GPU tensors.
    """
    return torch.sqrt(_as_float_tensor(x))


def _torch_erf(x: Any) -> Any:
    """Error function supporting real torch tensors and complex input."""
    if torch.is_tensor(x):
        if torch.is_complex(x):
            # torch.erf is not implemented for complex tensors.
            out = scipy_special.erf(x.detach().cpu().numpy())
            return torch.as_tensor(out, device=x.device, dtype=x.dtype)
        return torch.erf(_as_float_tensor(x))
    return scipy_special.erf(x)


def _torch_erfc(x: Any) -> Any:
    """Complementary error function with a complex-input fallback."""
    if torch.is_tensor(x):
        if torch.is_complex(x):
            out = scipy_special.erfc(x.detach().cpu().numpy())
            return torch.as_tensor(out, device=x.device, dtype=x.dtype)
        return torch.erfc(_as_float_tensor(x))
    return scipy_special.erfc(x)


def _torch_gamma(x: Any) -> Any:
    """Gamma function valid for negative and complex arguments.

    ``exp(lgamma(x))`` returns ``|Gamma(x)|`` and therefore has the wrong sign
    on ``(-1, 0)``, ``(-3, -2)``, ...  Negative arguments are handled with the
    reflection formula ``Gamma(x) = pi / (sin(pi x) Gamma(1 - x))``.
    """
    if not torch.is_tensor(x):
        return scipy_special.gamma(x)
    if torch.is_complex(x):
        out = scipy_special.gamma(x.detach().cpu().numpy())
        return torch.as_tensor(out, device=x.device, dtype=x.dtype)

    x_f = _as_float_tensor(x)
    positive = x_f > 0
    # Mask the arguments of both branches so that neither produces inf/nan
    # (which would otherwise poison the gradients of ``torch.where``).
    safe_pos = torch.where(positive, x_f, torch.ones_like(x_f))
    safe_neg = torch.where(positive, torch.ones_like(x_f), 1.0 - x_f)
    gamma_pos = torch.exp(torch.lgamma(safe_pos))
    gamma_neg = math.pi / (torch.sin(math.pi * x_f) * torch.exp(torch.lgamma(safe_neg)))
    return torch.where(positive, gamma_pos, gamma_neg)


def _torch_heaviside(x: Any, value_at_zero: Any = 0.0) -> torch.Tensor:
    """Heaviside step accepting SymPy's optional second argument."""
    x_f = _as_float_tensor(x)
    if torch.is_tensor(value_at_zero):
        values = value_at_zero.to(device=x_f.device, dtype=x_f.dtype)
        values = values.expand_as(x_f) if values.ndim else torch.full_like(x_f, float(values))
    else:
        values = torch.full_like(x_f, float(value_at_zero))
    return torch.heaviside(x_f, values)


@lru_cache(maxsize=1)
def _default_module_mapping() -> Dict[str, Any]:
    """SymPy → torch mapping used by ``lambdify`` (built once, then cached)."""
    return {
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
        "expm1": torch.expm1, "log1p": torch.log1p,
        # Special functions used by symbolic integral evaluation
        "erf": _torch_erf, "erfc": _torch_erfc, "gamma": _torch_gamma,
        "loggamma": torch.lgamma,
        # Roots / powers
        "sqrt": _torch_sqrt, "Pow": torch.pow,
        # Complex algebra
        "conjugate": torch.conj, "conj": torch.conj,
        "re": torch.real, "im": torch.imag, "arg": torch.angle,
        # Misc
        "Abs": torch.abs, "sign": torch.sign,
        "floor": torch.floor, "ceiling": torch.ceil,
        "Min": torch.minimum, "Max": torch.maximum,
        "Heaviside": _torch_heaviside,
        # Constants
        "pi": math.pi,
        "E": math.e,
    }


# ----------------------------------------------------------------------
# Quadrature rules (nodes + weights), cached per (a, b, N, device, dtype)
# ----------------------------------------------------------------------
@lru_cache(maxsize=64)
def _legendre_nodes_weights(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Gauss–Legendre nodes/weights on ``[-1, 1]`` (cached).

    ``numpy.polynomial.legendre.leggauss`` builds and diagonalises a companion
    matrix, which becomes the bottleneck for the large node counts that
    oscillatory (diffraction) integrands need: ~1.6 s at ``n=5001`` and ~105 s
    at ``n=20001``.  SciPy's Newton-based routine returns the same nodes and
    weights to ~1e-13 and is an order of magnitude faster there.
    """
    nodes, weights = scipy_special.roots_legendre(n)
    return np.asarray(nodes, dtype=np.float64), np.asarray(weights, dtype=np.float64)


@lru_cache(maxsize=256)
def _cached_rule(
    name: str, a: float, b: float, n: int, device: torch.device, dtype: torch.dtype
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Cached ``(coords, weights)`` for a built-in rule.

    The returned tensors are shared between calls and must be treated as
    read-only.
    """
    if name == "simpson":
        if n < 3 or n % 2 == 0:
            raise ValueError(f"Simpson rule requires odd N>=3; got N={n}")
        coords = torch.linspace(a, b, n, device=device, dtype=dtype)
        weights = torch.ones(n, device=device, dtype=dtype)
        weights[1:-1:2] = 4.0
        weights[2:-1:2] = 2.0
        weights *= (b - a) / (n - 1) / 3.0
        return coords, weights

    if name == "trapezoid":
        if n < 2:
            raise ValueError(f"Trapezoid rule requires N>=2; got N={n}")
        coords = torch.linspace(a, b, n, device=device, dtype=dtype)
        weights = torch.ones(n, device=device, dtype=dtype)
        weights[0] = 0.5
        weights[-1] = 0.5
        weights *= (b - a) / (n - 1)
        return coords, weights

    if name == "gauss-legendre":
        if n < 1:
            raise ValueError(f"Gauss-Legendre rule requires N>=1; got N={n}")
        nodes, weights = _legendre_nodes_weights(n)
        nodes_t = torch.as_tensor(nodes, device=device, dtype=dtype)
        weights_t = torch.as_tensor(weights, device=device, dtype=dtype)
        # Affine map from [-1, 1] to [a, b].
        coords = 0.5 * (b - a) * nodes_t + 0.5 * (b + a)
        return coords, 0.5 * (b - a) * weights_t

    raise ValueError(f"Unknown quadrature rule: {name}")  # pragma: no cover


#: Tensor-product grids are expensive to rebuild (``cartesian_prod`` + ``kron``),
#: so a few of them are kept around keyed by the rule/domain/N/device/dtype.
_GRID_CACHE: "OrderedDict[Any, Tuple[int, Callable]]" = OrderedDict()
_GRID_CACHE_MAXSIZE = 4
_GRID_CACHE_MAX_ELEMENTS = 1 << 23

#: Upper bound on ``n_quadrature_points * n_parameter_values`` per integrand
#: evaluation (2**27 float64 numbers ~ 1 GiB).  It is the budget used
#: internally by the separable / shift-invariant evaluators, and the budget
#: that ``chunk_size_params="auto"`` compares against when it decides whether
#: the vectorized solver should loop over chunks of the parameter mesh.  The
#: vectorized solver never chunks on its own: without an explicit
#: ``chunk_size_params`` it evaluates every parameter value in one call.
MAX_ELEMENTS_PER_EVALUATION = 1 << 27


_RULE_ALIASES: Dict[str, str] = {
    "simpson": "simpson",
    "trapezoid": "trapezoid",
    "trapezoidal": "trapezoid",
    "gauss-legendre": "gauss-legendre",
    "gauss_legendre": "gauss-legendre",
    "legendre": "gauss-legendre",
}


def clear_caches() -> None:
    """Drop every internal cache (quadrature nodes, grids, compiled functions).

    Useful to release device memory, or between benchmark runs.
    """
    _cached_rule.cache_clear()
    _legendre_nodes_weights.cache_clear()
    _default_module_mapping.cache_clear()
    _GRID_CACHE.clear()
    TorchSymPy._lambdify_cache.clear()


def _is_out_of_memory_error(exc: BaseException) -> bool:
    """Return ``True`` when ``exc`` reports a failed (CPU or GPU) allocation.

    Torch signals "cannot allocate" both as :class:`MemoryError` and as a plain
    :class:`RuntimeError` whose message mentions the allocator, so the message
    has to be inspected as well.
    """
    if isinstance(exc, MemoryError):
        return True
    cuda_oom = getattr(torch.cuda, "OutOfMemoryError", None)
    if cuda_oom is not None and isinstance(exc, cuda_oom):
        return True
    message = str(exc).lower()
    return (
        "out of memory" in message
        or "can't allocate memory" in message
        or "cannot allocate memory" in message
        or "defaultcpuallocator" in message
    )


def _resolve_device_dtype(
    device: Optional[torch.device], dtype: Optional[torch.dtype]
) -> Tuple[torch.device, torch.dtype, torch.dtype]:
    """Return ``(device, compute_dtype, param_dtype)``.

    ``compute_dtype`` is always a *real* dtype (quadrature nodes, weights and
    the real/imaginary outputs); ``param_dtype`` keeps the user's choice so
    that complex parameters survive.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(device)
    if dtype is None:
        dtype = torch.float32 if device.type == "cuda" else torch.float64
    return device, _real_dtype(dtype), dtype


# ----------------------------------------------------------------------
# Separable (oscillatory) multi-dimensional integrals
# ----------------------------------------------------------------------
#
# WHY THIS EXISTS
# ---------------
# A tensor-product quadrature spends its node budget ``N`` over *all* the
# integration axes: a 2D integral with ``N = 5001`` really uses only
# ``5001**(1/2) ~ 70`` nodes per axis.  That is fatal for diffraction
# integrands such as the Fresnel kernel
#
#     exp(i k/(2 z) ((x - x0)**2 + (y - y0)**2)),
#
# whose local frequency in the integration variable is ``|x - x0| / (lambda z)``:
# the number of fringes across the aperture grows as the screen distance ``z``
# shrinks, so the fixed 70 nodes per axis alias the oscillation - first at the
# large screen coordinates (a spurious "secondary pattern" along the four sides
# of the image) and, for small ``z``, everywhere.
#
# Raising ``N`` is not a cure, because a 2D rule needs the *square* of the nodes
# per axis (~1.4e7 points for the 40 mm case below).
#
# HOW IT IS FIXED HERE
# --------------------
# Almost every such integrand is *separable*: after expanding the phase, no
# single term mixes two integration variables, so
#
#     Int Int A(x0) B(y0) exp(i (F0 + Fx(x0) + Fy(y0))) dx0 dy0
#       = exp(i F0) [Int A exp(i Fx) dx0] [Int B exp(i Fy) dy0],
#
# and ``Int Int A B cos(F)`` / ``Int Int A B sin(F)`` are the real / imaginary
# part of that product (for real ``A``, ``B``).  The d-dimensional rule is
# replaced by ``d`` one-dimensional rules, so the *full* node budget is
# available on each axis, and the cost grows linearly - not as the d-th power -
# with the resolution.  Each 1D rule is then refined by doubling until the
# result stops changing, which removes the aliasing automatically whatever the
# distance ``z`` is.
#
# Two further details matter for the accuracy:
#   * the 1D factors depend on fewer parameters than the full screen mesh (the
#     x-integral only on x), so they are evaluated on the *unique* parameter
#     values and broadcast back - a 100x100 screen costs 100 evaluations, which
#     is what makes the refinement loop affordable;
#   * everything is computed in the real dtype resolved by
#     ``_resolve_device_dtype`` (float64 on the CPU); float32 caps these
#     oscillatory integrals at ~1e-3 relative accuracy.
#
# The decomposition is attempted for every ``Integral`` passed to
# :meth:`TorchSymPy.eval_numeric`; when it does not apply the original
# tensor-product solvers are used unchanged.

#: Nodes per axis used by the first pass of the separable evaluator (composite
#: Simpson needs an odd node count; ``2**m + 1`` keeps it odd while doubling).
SEPARABLE_START_NODES = 2 ** 11 + 1

#: Hard ceiling for the refinement loop (~2e6 nodes per axis).
SEPARABLE_MAX_NODES = 2 ** 21 + 1

#: Relative tolerance of the refinement loop.
SEPARABLE_TOL = 1e-9


def _is_finite_bound(bound) -> bool:
    """``True`` when an integration bound is a concrete finite number."""
    if bound in (oo, -oo):
        return False
    try:
        return math.isfinite(float(bound))
    except (TypeError, ValueError):
        return False


def _is_real_valued(expr) -> bool:
    """Conservative test: ``True`` when ``expr`` is real for real parameters."""
    return not expr.has(sympy_I, sympy_re, sympy_im, sympy_conjugate)


def _split_additive_by_vars(expr, variables):
    """Split ``expr`` into a variable-free part plus one part per variable.

    Returns ``(constant_part, {variable: part})``, or ``None`` when a term of
    the expanded expression involves two integration variables at once (the
    integrand is then genuinely coupled and cannot be separated).
    """
    parts: Dict[Symbol, List[Any]] = {variable: [] for variable in variables}
    constant: List[Any] = []
    for term in Add.make_args(sympy_expand(expr)):
        involved = [variable for variable in variables if term.has(variable)]
        if not involved:
            constant.append(term)
        elif len(involved) == 1:
            parts[involved[0]].append(term)
        else:
            return None
    return Add(*constant), {variable: Add(*terms) for variable, terms in parts.items()}


@dataclass
class _SeparableFactor:
    """One 1D factor ``Int amplitude * exp(i phase) dvariable``."""

    variable: Symbol
    lower: float
    upper: float
    amplitude: Any
    phase: Any


@dataclass
class _SeparablePlan:
    """Factorised form of a multi-dimensional integral.

    The value of the integral is::

        take(outer, take(trig, coeff * exp(i*const_phase) * prod(factors)))

    with ``trig`` in ``{None, "cos", "sin"}`` (real / imaginary part of the
    product) and ``outer`` in ``{None, "re", "im"}`` for an integrand that was
    wrapped in SymPy's ``re()`` / ``im()``.
    """

    coeff: Any
    const_phase: Any
    factors: List[_SeparableFactor]
    trig: Optional[str]
    outer: Optional[str]


def _separable_plan(integral: Integral) -> Optional[_SeparablePlan]:
    """Try to factorise ``integral`` into independent 1D integrals.

    ``None`` is returned whenever the integrand does not have the required
    shape, in which case the caller falls back to the tensor-product solvers.
    Recognised integrands are products of

    * factors depending on at most one integration variable,
    * ``exp(...)`` with an additively separable argument,
    * at most one ``cos(...)`` / ``sin(...)`` with an additively separable
      argument (its cofactors must then be real-valued),

    optionally wrapped in ``re(...)`` / ``im(...)``.
    """
    limits = list(integral.limits)
    if not limits:
        return None

    variables: List[Symbol] = []
    bounds: Dict[Symbol, Tuple[float, float]] = {}
    for limit in limits:
        if len(limit) != 3:
            return None
        variable, lower, upper = limit
        if not (_is_finite_bound(lower) and _is_finite_bound(upper)):
            return None
        variables.append(variable)
        bounds[variable] = (float(lower), float(upper))

    function = integral.function
    if function.has(Integral):
        return None

    outer: Optional[str] = None
    if isinstance(function, sympy_re):
        outer, function = "re", function.args[0]
    elif isinstance(function, sympy_im):
        outer, function = "im", function.args[0]

    amplitudes: Dict[Symbol, List[Any]] = {variable: [] for variable in variables}
    phases: Dict[Symbol, Any] = {variable: Integer(0) for variable in variables}
    coeff_factors: List[Any] = []
    const_phase: Any = Integer(0)
    trig: Optional[str] = None

    for factor in Mul.make_args(function):
        involved = [variable for variable in variables if factor.has(variable)]
        if not involved:
            coeff_factors.append(factor)
            continue

        if isinstance(factor, (sympy_cos, sympy_sin)):
            if trig is not None:
                return None
            split = _split_additive_by_vars(factor.args[0], variables)
            if split is None:
                return None
            const_phase, phases = split
            if not all(_is_real_valued(part) for part in (const_phase, *phases.values())):
                return None
            trig = "cos" if isinstance(factor, sympy_cos) else "sin"
            continue

        if isinstance(factor, sympy_exp):
            split = _split_additive_by_vars(factor.args[0], variables)
            if split is None:
                return None
            constant_part, variable_parts = split
            if constant_part != 0:
                coeff_factors.append(sympy_exp(constant_part))
            for variable, part in variable_parts.items():
                if part != 0:
                    amplitudes[variable].append(sympy_exp(part))
            continue

        if len(involved) != 1:
            return None
        amplitudes[involved[0]].append(factor)

    factors = [
        _SeparableFactor(
            variable=variable,
            lower=bounds[variable][0],
            upper=bounds[variable][1],
            amplitude=Mul(*amplitudes[variable]) if amplitudes[variable] else Integer(1),
            phase=phases[variable],
        )
        for variable in variables
    ]
    coeff = Mul(*coeff_factors) if coeff_factors else Integer(1)

    if trig is not None:
        # Re/Im may only be pulled out of the product when everything else is
        # real; otherwise the factorisation would silently change the result.
        cofactors = [coeff] + [factor.amplitude for factor in factors]
        if not all(_is_real_valued(cofactor) for cofactor in cofactors):
            return None

    return _SeparablePlan(
        coeff=coeff, const_phase=const_phase, factors=factors, trig=trig, outer=outer
    )


# ----------------------------------------------------------------------
# Shift-invariant (convolution) integrands: Rayleigh-Sommerfeld & friends
# ----------------------------------------------------------------------
#
# The separable decomposition above cannot help the Rayleigh-Sommerfeld
# integrand
#
#     E(x,y) = 1/(i lambda) Int Int z exp(i k r)/r^2 dx0 dy0,
#     r = sqrt((x-x0)^2 + (y-y0)^2 + z^2),
#
# because the square root couples x0 and y0: no rewriting turns it into a
# product of a function of x0 and a function of y0.  It suffers from exactly
# the same aliasing as the Fresnel one, though - the phase k r advances by
# ~|x-x0|/(lambda z) cycles per unit length, i.e. 47 cycles across the aperture
# at z = 400 mm and 475 at z = 40 mm for the Abedin2005Fig6a configuration -
# and a fixed budget of N = 5001 tensor-product nodes (~70 per axis) samples
# far too coarsely.  Because the local frequency grows with |x-x0| the aliasing
# appears first at the largest screen coordinates, i.e. as a spurious pattern
# along the four sides of the image, and takes over completely at short z.
#
# What this kernel *does* have is shift invariance: the integrand depends on
# the integration variables and the screen coordinates only through the
# differences u = x - x0 and v = y - y0.  Substituting them,
#
#     Int_{a1}^{a2} Int_{b1}^{b2} h(x-x0, y-y0) dy0 dx0
#         = Int_{x-a2}^{x-a1} Int_{y-b2}^{y-b1} h(u,v) dv du
#         = Phi(x-a1, y-b1) - Phi(x-a2, y-b1)
#           - Phi(x-a1, y-b2) + Phi(x-a2, y-b2),
#     Phi(a,b) = Int^{a} Int^{b} h(u,v) dv du,
#
# so *every* screen point is a difference of four values of one cumulative
# integral of the kernel.  Splitting the u-axis at all the distinct corner
# abscissae x-a1, x-a2 (200 of them for a 100-point screen) and integrating
# each panel with Gauss-Legendre gives all of those values at once: the panel
# integrals are computed on a (panels x nodes)^2 grid and turned into Phi by a
# 2D cumulative sum.  The kernel is therefore evaluated O(n_nodes^2) times in
# total instead of O(n_nodes^2) times *per screen point*, which is what makes a
# properly resolved rule affordable; the node count per panel is doubled until
# two successive results agree, so the resolution again adapts to lambda, z and
# the screen size instead of being fixed.
#
# The plan below recognises this structure symbolically (for one or two
# integration variables) by substituting x0 -> x - u and checking that the
# screen coordinate disappears from the integrand.  Anything else falls back to
# the original tensor-product solvers.

#: Gauss-Legendre nodes per panel used by the first pass of the shift-invariant
#: evaluator.
SHIFT_INVARIANT_START_NODES = 8

#: Hard ceiling for its refinement loop (nodes per panel).
SHIFT_INVARIANT_MAX_NODES = 1024

#: Relative tolerance of its refinement loop.
SHIFT_INVARIANT_TOL = 1e-9

#: Largest integrand (in ``count_ops``) on which the recognition below is
#: allowed to call ``simplify``.
SHIFT_INVARIANT_MAX_OPS = 200

#: Corner abscissae closer than this (relative to their magnitude) are treated
#: as one panel boundary, so that round-off in ``x - a`` does not create a
#: swarm of degenerate panels.  It is a few units in the last place: merging
#: needs no more than that, while moving a boundary by ``tol`` moves the
#: result by ``tol`` times the size of the integrand.
SHIFT_INVARIANT_BREAKPOINT_TOL = 1e-15


@dataclass
class _ShiftAxis:
    """One integration axis of a shift-invariant integrand.

    The integration variable enters the integrand only through
    ``coordinate = sign * (parameter - variable)``; the panel boundaries needed
    for that axis are ``sign * parameter + offset_lower`` and
    ``sign * parameter + offset_upper``.
    """

    variable: Symbol
    parameter: Symbol
    coordinate: Symbol
    sign: int
    offset_lower: float
    offset_upper: float


@dataclass
class _ShiftInvariantPlan:
    """Convolution form of an integral: kernel plus one axis description."""

    kernel: Any
    axes: List[_ShiftAxis]


def _shift_invariant_plan(integral: Integral, parameters) -> Optional[_ShiftInvariantPlan]:
    """Recognise ``Int h(p1 - v1, p2 - v2) dv1 dv2`` (a convolution).

    ``parameters`` are the free symbols of the surrounding expression that may
    play the role of the screen coordinates.  ``None`` is returned when the
    integrand does not have that shape, or when it has more than two
    integration variables, or when a bound is not a finite number - the caller
    then falls back to the tensor-product solvers.
    """
    limits = list(integral.limits)
    if not 1 <= len(limits) <= 2:
        return None

    variables: List[Symbol] = []
    bounds: List[Tuple[float, float]] = []
    for limit in limits:
        if len(limit) != 3:
            return None
        variable, lower, upper = limit
        if not (_is_finite_bound(lower) and _is_finite_bound(upper)):
            return None
        variables.append(variable)
        bounds.append((float(lower), float(upper)))

    function = integral.function
    if function.has(Integral):
        return None

    candidates = [
        symbol
        for symbol in parameters
        if function.has(symbol) and symbol not in variables
    ]
    if len(candidates) < len(variables) or len(candidates) > 4:
        return None

    for assignment in permutations(candidates, len(variables)):
        for signs in product((1, -1), repeat=len(variables)):
            substitution = {}
            coordinates = []
            for variable, parameter, sign in zip(variables, assignment, signs):
                coordinate = Dummy("u_" + variable.name, real=True)
                coordinates.append(coordinate)
                # sign=+1: u = parameter - variable; sign=-1: u = variable - parameter.
                substitution[variable] = (
                    parameter - coordinate if sign > 0 else parameter + coordinate
                )
            kernel = function.subs(substitution, simultaneous=True)
            if any(kernel.has(parameter) for parameter in assignment):
                # (x - x0)^2 is often stored expanded, e.g. inside a sqrt, so
                # the parameter only cancels after expanding; simplify() is the
                # last resort and is skipped for large expressions, where it
                # could cost more than the integration itself.
                kernel = sympy_expand(kernel)
            if any(kernel.has(parameter) for parameter in assignment):
                if count_ops(kernel) > SHIFT_INVARIANT_MAX_OPS:
                    continue
                kernel = sympy_simplify(kernel)
                if any(kernel.has(parameter) for parameter in assignment):
                    continue

            axes = []
            for variable, parameter, sign, coordinate, (lower, upper) in zip(
                variables, assignment, signs, coordinates, bounds
            ):
                # sign=+1: u runs over [parameter-upper, parameter-lower];
                # sign=-1: u runs over [lower-parameter, upper-parameter].
                offsets = (-upper, -lower) if sign > 0 else (lower, upper)
                axes.append(
                    _ShiftAxis(
                        variable=variable,
                        parameter=parameter,
                        coordinate=coordinate,
                        sign=sign,
                        offset_lower=offsets[0],
                        offset_upper=offsets[1],
                    )
                )
            return _ShiftInvariantPlan(kernel=kernel, axes=axes)
    return None


def _panel_boundaries(corners: np.ndarray, tol: float = SHIFT_INVARIANT_BREAKPOINT_TOL):
    """Distinct panel boundaries of ``corners`` plus the index of each corner."""
    scale = max(1.0, float(np.abs(corners).max()) if corners.size else 1.0)
    step = tol * scale
    rounded = np.round(corners / step) * step
    boundaries, inverse = np.unique(rounded, return_inverse=True)
    return boundaries, inverse


def _panel_nodes_weights(boundaries: torch.Tensor, n_nodes: int):
    """Gauss-Legendre nodes/weights of every panel: two ``(panels, n_nodes)`` tensors."""
    reference_nodes, reference_weights = _legendre_nodes_weights(n_nodes)
    nodes = torch.as_tensor(
        reference_nodes, device=boundaries.device, dtype=boundaries.dtype
    )
    weights = torch.as_tensor(
        reference_weights, device=boundaries.device, dtype=boundaries.dtype
    )
    lower, upper = boundaries[:-1], boundaries[1:]
    middle = 0.5 * (lower + upper)
    half = 0.5 * (upper - lower)
    return (
        middle.unsqueeze(1) + half.unsqueeze(1) * nodes.unsqueeze(0),
        half.unsqueeze(1) * weights.unsqueeze(0),
    )


def _broadcast_values(values, shape, *, device, dtype) -> torch.Tensor:
    """Tensor of shape ``shape`` from whatever an integrand callable returned."""
    if not torch.is_tensor(values):
        values = torch.as_tensor(values, device=device, dtype=dtype)
    return torch.broadcast_to(values, shape)


@dataclass
class TorchExpr:
    """Pre-built torch function together with its finite integration domain.

    This is the object returned by :meth:`TorchSymPy.torchify` when integration
    limits are supplied (directly or via a SymPy ``Integral``).  It stores both
    the numerical callable and all metadata needed by the integration helpers
    in this module.

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
        Finite box for the integration domain.  Each entry is ``[a, b]`` for
        one dimension after change-of-variables.  This is what the integrators
        use as their integration range.
    dim : int
        Number of integration variables (i.e. the dimensionality of the
        domain).
    n_params : int
        Number of additional symbolic parameters expected by ``func``.
    sympy_integrand : Any
        The SymPy integrand *after* all symbolic substitutions and Jacobians
        have been applied.  Kept for debugging / inspection.
    sympy_integral : Any
        The original SymPy ``Integral`` if there was one, else ``None``.
    variables : list[sympy.Symbol] | None
        Integration variables used by the numerical object after any
        change-of-variables has been applied.  For infinite limits these are
        the transformed symbols such as ``t_x`` rather than the original
        symbolic variable ``x``.
    """

    func: Callable                      # f(new_var_0, …, new_var_n, param_0, …, param_m)
    domain: List[List[float]]           # finite box for the quadrature
    dim: int                            # number of integration variables
    n_params: int                       # number of extra parameters
    sympy_integrand: Any = None         # reduced sympy integrand (if any)
    sympy_integral: Any = None          # original sympy integral (if any)
    variables: Optional[List[Symbol]] = None  # sympy variables after change of variables

    # ------------------------------------------------------------------
    # Quadrature rules
    # ------------------------------------------------------------------
    @staticmethod
    def _rule_simpson(a: float, b: float, N: int, *, device=None, dtype=None):
        """Composite Simpson nodes and weights on ``[a, b]`` (odd ``N >= 3``)."""
        return _cached_rule("simpson", float(a), float(b), int(N), device, dtype)

    @staticmethod
    def _rule_trapezoid(a: float, b: float, N: int, *, device=None, dtype=None):
        """Composite trapezoid nodes and weights on ``[a, b]`` (``N >= 2``)."""
        return _cached_rule("trapezoid", float(a), float(b), int(N), device, dtype)

    @staticmethod
    def _rule_gauss_legendre(a: float, b: float, N: int, *, device=None, dtype=None):
        """Gauss–Legendre nodes and weights mapped from ``[-1, 1]`` to ``[a, b]``."""
        return _cached_rule("gauss-legendre", float(a), float(b), int(N), device, dtype)

    @staticmethod
    def _resolve_quadrature_rule(method: Union[str, QuadratureRule, None]) -> QuadratureRule:
        """Map a rule name (or a user callable) to a ``(a, b, N) -> (x, w)`` callable."""
        if callable(method):
            return method

        method_name = "simpson" if method is None else str(method).strip().lower()
        rule_name = _RULE_ALIASES.get(method_name)
        if rule_name is None:
            raise ValueError(
                "Unknown integration method. Expected one of "
                "{'simpson', 'trapezoid', 'gauss-legendre'} or a callable. "
                f"Got: {method}"
            )

        def rule(a, b, N, *, device=None, dtype=None):
            return _cached_rule(rule_name, float(a), float(b), int(N), device, dtype)

        return rule

    # ------------------------------------------------------------------
    # Input normalisation helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _is_scalar_like(x: Any) -> bool:
        """Return ``True`` for scalar-like values.

        This is intentionally permissive: besides builtin numeric types, it
        treats SymPy/mpmath numeric scalars as scalar-like if they can be
        coerced to ``float`` or ``complex``.
        """
        if torch.is_tensor(x):
            return x.ndim == 0
        if isinstance(x, (int, float, complex, np.number)):
            return True
        if isinstance(x, np.ndarray):
            return x.ndim == 0
        for cast in (float, complex):
            try:
                cast(x)
                return True
            except (TypeError, ValueError):
                continue
        return False

    @staticmethod
    def _to_tensor(value: Any, *, device=None, dtype=None) -> torch.Tensor:
        """Best-effort conversion of scalars / arrays / SymPy numbers to a tensor.

        ``torch.as_tensor`` cannot infer a dtype for some scalar-like types
        (``mpmath.mpf``, ``sympy.Float``, object-dtype numpy arrays), so those
        are coerced to builtin ``float``/``complex`` first.
        """
        if torch.is_tensor(value):
            return value.to(device=device, dtype=dtype) if dtype is not None else value.to(device=device)

        try:  # fast path
            return torch.as_tensor(value, device=device, dtype=dtype)
        except (TypeError, RuntimeError, ValueError):
            pass

        array = value if isinstance(value, np.ndarray) else None
        if array is None:
            try:
                array = np.asarray(value)
            except (TypeError, ValueError):
                array = None

        if isinstance(array, np.ndarray) and array.ndim > 0:
            if array.dtype == object:
                # Object arrays typically hold SymPy/mpmath numbers.
                try:
                    array = array.astype(np.float64)
                except (TypeError, ValueError):
                    array = np.vectorize(complex, otypes=[np.complex128])(array)
            return torch.as_tensor(array, device=device)

        # Scalar fallback: coerce to builtin numeric types (complex first so
        # that the imaginary part is never silently dropped).
        try:
            return torch.as_tensor(complex(value), device=device)
        except (TypeError, ValueError):
            return torch.as_tensor(float(value), device=device, dtype=dtype)

    @staticmethod
    def _normalize_params_values(
        params_values: Any, n_params: int, *, device, dtype
    ) -> Tuple[torch.Tensor, Tuple[int, ...]]:
        """Normalize parameter inputs into a 2D tensor and its batch shape.

        Accepted input formats are:

        - ``None`` when ``n_params == 0``,
        - a tensor of shape ``(..., n_params)``, or
        - a list/tuple of length ``n_params`` whose elements are scalars or
          tensors with the same shape.

        Returns
        -------
        params_flat : torch.Tensor
            Shape ``(batch_size, n_params)``.
        batch_shape : tuple[int, ...]
            Shape the results have to be reshaped to.
        """
        if params_values is None:
            if n_params != 0:
                raise ValueError(f"Expected {n_params} params; got None")
            return torch.empty((1, 0), device=device, dtype=dtype), ()

        if torch.is_tensor(params_values) or isinstance(params_values, np.ndarray):
            params_tensor = TorchExpr._to_tensor(params_values, device=device)
            if params_tensor.ndim == 0 or params_tensor.shape[-1] != n_params:
                raise ValueError(
                    "params_values last dimension must be n_params="
                    f"{n_params}; got shape={tuple(params_tensor.shape)}"
                )
            batch_shape = tuple(params_tensor.shape[:-1])
            batch_size = math.prod(batch_shape) if batch_shape else 1
            params_flat = params_tensor.reshape(batch_size, n_params).to(device=device, dtype=dtype)
            return params_flat, batch_shape

        if not isinstance(params_values, (list, tuple)):
            params_values = [params_values]
        if len(params_values) != n_params:
            raise ValueError(f"Expected {n_params} params; got {len(params_values)}")

        tensors = [TorchExpr._to_tensor(value, device=device) for value in params_values]

        batch_shape: Tuple[int, ...] = ()
        for tensor in tensors:
            if tensor.ndim > 0:
                batch_shape = tuple(tensor.shape)
                break

        aligned = []
        for tensor in tensors:
            if tensor.ndim == 0:
                aligned.append(tensor.expand(batch_shape) if batch_shape else tensor)
            elif tuple(tensor.shape) != batch_shape:
                raise ValueError(
                    "All parameter tensors must have the same shape; got "
                    f"{tuple(tensor.shape)} vs {batch_shape}"
                )
            else:
                aligned.append(tensor)

        params_tensor = torch.stack(aligned, dim=-1) if aligned else torch.empty(
            batch_shape + (0,), device=device, dtype=dtype
        )
        if dtype is not None:
            params_tensor = params_tensor.to(dtype=dtype)
        batch_size = math.prod(batch_shape) if batch_shape else 1
        return params_tensor.reshape(batch_size, n_params), batch_shape

    @staticmethod
    def _align_values(values: Any, n_points: int, n_batch: int, *, device) -> torch.Tensor:
        """Coerce an integrand result to the shape ``(n_points, n_batch)``.

        The integrand is called with points of shape ``(n_points, 1)`` and
        parameters of shape ``(1, n_batch)``, so a correct result broadcasts to
        ``(n_points, n_batch)``.  Constant or partially-broadcast results (for
        instance an integrand that does not depend on the parameters) are
        expanded here instead of being rejected.
        """
        values = values if torch.is_tensor(values) else TorchExpr._to_tensor(values, device=device)
        if values.device != device:
            values = values.to(device=device)

        target = (n_points, n_batch)
        if values.ndim == 0:
            values = values.reshape(1, 1)
        elif values.ndim == 1:
            if values.shape[0] == n_points:
                values = values.unsqueeze(1)
            elif values.shape[0] == n_batch:
                values = values.unsqueeze(0)
            elif values.shape[0] == 1:
                values = values.reshape(1, 1)
            else:
                raise RuntimeError(
                    f"Unexpected integrand output shape {tuple(values.shape)}; expected {target}"
                )
        elif values.ndim == 2:
            rows, cols = values.shape
            row_ok = rows in (n_points, 1)
            col_ok = cols in (n_batch, 1)
            if not (row_ok and col_ok) and (rows, cols) == (n_batch, n_points):
                # Only transpose when the result cannot be read as
                # (points, batch); doing it unconditionally corrupted the
                # results whenever n_points == n_batch.
                values = values.t()
        else:
            raise RuntimeError(
                f"Unexpected integrand output shape {tuple(values.shape)}; expected {target}"
            )

        if tuple(values.shape) != target:
            try:
                values = values.expand(target).contiguous()
            except RuntimeError as exc:
                raise RuntimeError(
                    f"Unexpected integrand output shape {tuple(values.shape)}; expected {target}"
                ) from exc
        return values

    @staticmethod
    def _make_grid_chunker(
        coords_1d: Sequence[torch.Tensor],
        weights_1d: Sequence[torch.Tensor],
        *,
        materialize_limit: int = 1 << 23,
    ) -> Tuple[int, Callable[[int, int], Tuple[torch.Tensor, torch.Tensor]]]:
        """Build a ``(start, stop) -> (points, weights)`` accessor for the grid.

        Small tensor-product grids are materialised once (fastest); larger ones
        are generated on the fly by index arithmetic, so a ``dim``-dimensional
        rule with ``N`` points per axis never needs ``N**dim * dim`` numbers in
        memory at the same time.
        """
        dim = len(coords_1d)
        sizes = [int(c.numel()) for c in coords_1d]
        n_points = math.prod(sizes)
        device = coords_1d[0].device

        if n_points * dim <= materialize_limit:
            if dim == 1:
                grid = coords_1d[0].unsqueeze(1)
                weights = weights_1d[0]
            else:
                grid = torch.cartesian_prod(*coords_1d)
                weights = weights_1d[0]
                for weight_1d in weights_1d[1:]:
                    weights = torch.kron(weights, weight_1d)

            def take(start: int, stop: int):
                return grid[start:stop], weights[start:stop]

            return n_points, take

        # ``cartesian_prod`` varies the last axis fastest; mirror that here.
        strides = [1] * dim
        for axis in range(dim - 2, -1, -1):
            strides[axis] = strides[axis + 1] * sizes[axis + 1]

        def take_lazy(start: int, stop: int):
            flat_index = torch.arange(start, stop, device=device)
            columns = []
            weights = None
            for axis in range(dim):
                index = torch.div(flat_index, strides[axis], rounding_mode="floor") % sizes[axis]
                columns.append(coords_1d[axis][index])
                weight = weights_1d[axis][index]
                weights = weight if weights is None else weights * weight
            return torch.stack(columns, dim=1), weights

        return n_points, take_lazy

    @staticmethod
    def _cached_grid_chunker(key, coords_1d, weights_1d):
        """``_make_grid_chunker`` with a small cache of materialised grids."""
        cached = _GRID_CACHE.get(key)
        if cached is not None:
            _GRID_CACHE.move_to_end(key)
            return cached

        n_points, take = TorchExpr._make_grid_chunker(coords_1d, weights_1d)
        if n_points * len(coords_1d) <= _GRID_CACHE_MAX_ELEMENTS:
            _GRID_CACHE[key] = (n_points, take)
            while len(_GRID_CACHE) > _GRID_CACHE_MAXSIZE:
                _GRID_CACHE.popitem(last=False)
        return n_points, take

    # ------------------------------------------------------------------
    # Convenience methods: numeric integration on this TorchExpr
    # ------------------------------------------------------------------
    def torchquad_integrate(self, params_values=None, method=None, N: int = 21, dtype=None):
        """Integrate this ``TorchExpr`` using the torchquad-based integrator.

        This is the simplest integration entry point.  It is intended for
        low-dimensional problems where a single torchquad integral is
        sufficient.

        Parameters
        ----------
        params_values : list | tuple | torch.Tensor | None, optional
            Numerical values for the symbolic parameters:

            - ``None`` if there are no parameters,
            - a Python sequence (list/tuple) of numbers/tensors, or
            - a tensor with one entry per parameter.
        method : torchquad integrator instance or None, optional
            If ``None``, a default ``Simpson`` integrator from torchquad is
            constructed internally.
        N : int, optional
            Resolution parameter passed through to torchquad.  Larger values
            increase accuracy but also the number of function evaluations.
        dtype : torch.dtype or None, optional
            Floating dtype used for internal computations.  This controls how
            domain points and parameter values are coerced before evaluation.
            If ``None``, torchquad's defaults are used.

        Returns
        -------
        re, im : torch.Tensor
            Real and imaginary parts of the integral.  If the integrand is
            real-valued, ``im`` is a zero tensor.

        USAGE
        =====
        1. Simple 1D integral without parameters::

            ````python
            import torchsympy
            from sympy import symbols, Integral, exp, oo

            x = symbols("x", real=True)
            integral_expr = Integral(exp(-x**2), (x, -oo, oo))

            lt = torchsympy.TorchSymPy()
            texpr = lt.torchify(integral_expr)
            re, im = texpr.torchquad_integrate(N=121)
            ````

        2. Integral with parameters (e.g. Fourier transform)::

            ````python
            import torchsympy
            from sympy import symbols, Integral, exp, I, oo

            x, k = symbols("x k", real=True)
            integral_expr = Integral(exp(-x**2) * exp(I * k * x), (x, -oo, oo))

            lt = torchsympy.TorchSymPy()
            texpr = lt.torchify(integral_expr)

            # Evaluate at k = 0.5
            re, im = texpr.torchquad_integrate(params_values=[0.5], N=121)
            ````
        """
        from torchquad import Simpson

        if method is None:
            method = Simpson()

        # The ``is None`` test is deliberate: the truthiness of a tensor with
        # more than one element raises.
        if params_values is None:
            param_values: List[Any] = []
        elif torch.is_tensor(params_values) or isinstance(params_values, np.ndarray):
            param_values = list(params_values.reshape(-1)) if params_values.ndim else [params_values]
        elif isinstance(params_values, (list, tuple)):
            param_values = list(params_values)
        else:
            param_values = [params_values]

        if len(param_values) != self.n_params:
            raise ValueError(f"Expected {self.n_params} params; got {len(param_values)}")

        # Convert the parameters once, not on every integrand evaluation.
        converted: Dict[torch.device, List[torch.Tensor]] = {}

        def params_on(device: torch.device) -> List[torch.Tensor]:
            cached = converted.get(device)
            if cached is None:
                cached = [
                    self._to_tensor(value, device=device, dtype=dtype) for value in param_values
                ]
                converted[device] = cached
            return cached

        def integrand(domain_points: torch.Tensor) -> torch.Tensor:
            points = domain_points.to(dtype=dtype) if dtype is not None else domain_points
            values = self.func(
                *(points[:, i] for i in range(points.shape[1])), *params_on(points.device)
            )
            out = values if torch.is_tensor(values) else torch.as_tensor(values, device=points.device)
            if dtype is not None and out.is_floating_point():
                out = out.to(dtype=dtype)
            return out

        result = method.integrate(
            integrand,
            dim=self.dim,
            N=N,
            integration_domain=self.domain,
        )
        if torch.is_complex(result):
            return result.real, result.imag
        return result, torch.zeros_like(result)

    def torchquad_integrate_vectorized(
        self,
        params_values=None,
        method=None,
        N: int = 21,
        chunk_size_params: Union[int, str, None] = None,
    ):
        """Vectorized integrator that passes raw parameter tensors to the integrand.

        This method skips nested Python loops and delegates all broadcasting to
        your integrand function.  The integrand receives the grid points of the
        integration variables (reshaped so that they occupy the leading axis)
        along with whatever raw, potentially multi-dimensional parameter
        tensors you pass in ``params_values``.

        It is your integrand's responsibility to combine these shapes
        mathematically without a size mismatch, either by reshaping manually
        (e.g. ``x.view(-1, 1, 1)``) or by using ``torch.einsum``.

        Parameters
        ----------
        params_values : iterable of torch.Tensor or None
            Raw parameter tensors passed directly to the integrand.
        method : torchquad integrator instance or None
            Integration method (e.g. ``Simpson()``, ``MonteCarlo()``,
            ``Boole()``).  If ``None``, defaults to ``Simpson()``.
        N : int
            Number of integration points.  Note that torchquad spreads ``N``
            over all ``dim`` axes, so a 2D integral with ``N=5001`` really uses
            only about ``5001 ** (1 / 2) ~ 70`` nodes per axis.
        chunk_size_params : int, ``"auto"`` or None, optional
            ``None`` (the default) keeps this method as simple as it has always
            been: a single vectorized call with every parameter value evaluated
            at once, and no Python loop.

            Pass a positive integer to evaluate at most that many parameter
            values at a time.  The integrand is evaluated on a tensor of shape
            ``(n_quadrature_points, *param_shape)``, so a very large parameter
            mesh (e.g. a 1000x1000 screen) may not fit in memory in one go;
            chunking caps the peak usage at
            ``n_quadrature_points * chunk_size_params`` elements, at the price
            of a Python loop over the chunks.

            Pass the string ``"auto"`` to have the chunk size picked
            automatically: everything at once while
            ``n_quadrature_points * n_parameter_values`` stays below
            :data:`MAX_ELEMENTS_PER_EVALUATION`, chunked above it.  The result
            is the same in all three cases; only peak memory and speed differ.

        Returns
        -------
        torch.Tensor
            The integration result.  For scalar-valued integrands the returned
            tensor is shape-normalized to match the batched API output shape.

        USAGE
        =====
        1. Broadcasting with manual reshaping (``.view()`` / ``.unsqueeze()``)::

            ````python
            import torchsympy
            lt = torchsympy.TorchSymPy()

            def integrand(x, A_grid, B_grid):
                # x is (101,) from torchquad, A_grid/B_grid are (10, 4000)
                x = x.view(-1, 1, 1)        # (101, 1, 1)
                A = A_grid.unsqueeze(0)     # (1, 10, 4000)
                B = B_grid.unsqueeze(0)     # (1, 10, 4000)
                return torch.sqrt(x * A * B)  # (101, 10, 4000)

            texpr = lt.torchify_callable(integrand, domain=[[0.0, 1.0]], n_params=1)
            texpr.torchquad_integrate_vectorized(params_values=[A_grid, B_grid])
            ````

        2. Broadcasting natively with ``torch.einsum``::

            ````python
            def integrand(x, grid_param):
                # x is (101,), grid_param is (10, 4000) → output (101, 10, 4000)
                return torch.sin(torch.einsum("i,jk->ijk", x, grid_param))

            texpr = lt.torchify_callable(integrand, domain=[[0.0, 1.0]], n_params=1)
            texpr.torchquad_integrate_vectorized(params_values=[grid_param])
            ````
        """
        if method is None:
            from torchquad import Simpson

            method = Simpson()

        params = list(params_values) if params_values is not None else []
        chunk_size_params = self._resolve_chunk_size_params(chunk_size_params, params, N)

        def integrate_with(local_params: Sequence[Any]):
            def integrand(domain_points: torch.Tensor):
                device = domain_points.device
                tensors = [self._to_tensor(p, device=device) for p in local_params]

                # Add one trailing singleton per parameter dimension so that the
                # quadrature points always sit on the leading axis.
                max_param_ndim = max((p.ndim for p in tensors), default=0)
                target_shape = [-1] + [1] * max_param_ndim
                var_args = [domain_points[:, i].view(*target_shape) for i in range(self.dim)]

                return self.func(*var_args, *tensors)

            return method.integrate(
                integrand, dim=self.dim, N=N, integration_domain=self.domain
            )

        # No chunking unless the caller asked for it: with chunk_size_params
        # left at None this is one plain vectorized call, exactly as before the
        # chunking option existed.
        if chunk_size_params is None:
            result = integrate_with(params)
        else:
            result = self._integrate_in_param_chunks(params, integrate_with, chunk_size_params)

        if torch.is_tensor(result) and result.ndim > 0 and result.shape[-1] == 1:
            result = result.squeeze(-1)
        return result

    @staticmethod
    def _is_auto_chunk_size(chunk_size_params: Any) -> bool:
        """True for the literal string ``"auto"`` (case/space insensitive)."""
        return (
            isinstance(chunk_size_params, str)
            and chunk_size_params.strip().lower() == "auto"
        )

    def _resolve_chunk_size_params(
        self, chunk_size_params: Any, params: Sequence[Any], N: int
    ) -> Optional[int]:
        """Turn the user's ``chunk_size_params`` into a chunk size or ``None``.

        ``None`` means "do not chunk", i.e. the plain single-call behaviour of
        the vectorized solver; a positive integer is used as given; only the
        literal string ``"auto"`` asks for a chunk size to be chosen here.
        """
        if chunk_size_params is None:
            return None
        if self._is_auto_chunk_size(chunk_size_params):
            return self._auto_chunk_size_params(params, N)
        if isinstance(chunk_size_params, bool) or not isinstance(
            chunk_size_params, (int, np.integer)
        ):
            raise ValueError(
                "chunk_size_params must be a positive integer, 'auto' or None; "
                f"got {chunk_size_params!r}"
            )
        if chunk_size_params <= 0:
            raise ValueError(f"chunk_size_params must be positive; got {chunk_size_params}")
        return int(chunk_size_params)

    def _auto_chunk_size_params(self, params: Sequence[Any], N: int) -> Optional[int]:
        """Chunk size for ``chunk_size_params="auto"``, or ``None`` for "no chunking".

        The integrand is evaluated on ``n_quadrature_points x n_parameter_values``
        numbers.  Rather than letting a large parameter mesh fail with an
        allocation error, cap that product at
        :data:`MAX_ELEMENTS_PER_EVALUATION`; the result is unchanged, only the
        peak memory (and a Python loop) differ.
        """
        batch_sizes = [
            int(np.prod(np.shape(p)))
            for p in params
            if not self._is_scalar_like(p) and np.ndim(p) > 0
        ]
        if not batch_sizes:
            return None
        batch_size = max(batch_sizes)

        # torchquad spreads N over the dim axes; this only needs to be a rough
        # estimate, since the chunk size affects memory and not the result.
        points = max(1, int(math.ceil(N ** (1.0 / max(self.dim, 1)))) ** self.dim)
        if points * batch_size <= MAX_ELEMENTS_PER_EVALUATION:
            return None

        chunk = max(1, MAX_ELEMENTS_PER_EVALUATION // points)
        logger.debug(
            "vectorized quadrature: {} parameter values x ~{} quadrature points "
            "exceeds the {} element budget; chunking into groups of {}",
            batch_size, points, MAX_ELEMENTS_PER_EVALUATION, chunk,
        )
        return chunk

    def _integrate_in_param_chunks(
        self,
        params: Sequence[Any],
        integrate_with: Callable[[Sequence[Any]], Any],
        chunk_size_params: int,
    ) -> torch.Tensor:
        """Apply ``integrate_with`` to slices of the parameter mesh, then stitch.

        All non-scalar parameters must share one common shape; they are
        flattened, sliced into chunks of at most ``chunk_size_params`` values,
        and the per-chunk results are concatenated and reshaped back to that
        common shape.
        """
        tensors = [self._to_tensor(p) for p in params]
        batch_shapes = {tuple(t.shape) for t in tensors if t.ndim > 0}
        if not batch_shapes:
            return integrate_with(tensors)
        if len(batch_shapes) > 1:
            raise ValueError(
                "chunk_size_params requires all non-scalar parameters to have the "
                f"same shape; got {sorted(batch_shapes)}"
            )

        batch_shape = batch_shapes.pop()
        batch_size = math.prod(batch_shape)
        flat = [t.reshape(-1) if t.ndim > 0 else t for t in tensors]

        pieces: List[torch.Tensor] = []
        for start in range(0, batch_size, chunk_size_params):
            stop = min(batch_size, start + chunk_size_params)
            chunk_params = [t[start:stop] if t.ndim > 0 else t for t in flat]
            piece = integrate_with(chunk_params)
            piece = piece if torch.is_tensor(piece) else self._to_tensor(piece)
            pieces.append(piece.reshape(-1))

        return torch.cat(pieces).reshape(batch_shape)

    def torch_integrate_batched(
        self,
        *,
        params_values=None,
        method: Union[str, QuadratureRule] = "simpson",
        N: int = 121,
        chunk_size_params: int = 256,
        chunk_size_points: Optional[int] = None,
        device=None,
        dtype=None,
    ):
        """Batched tensor-product quadrature integration for this ``TorchExpr``.

        This method is intended for situations where you want to evaluate the
        same integral for many different parameter values (for example, on a 2D
        or 3D mesh).

        Parameters
        ----------
        params_values : tensor | list | tuple | None
            Numerical parameters; see :meth:`_normalize_params_values` for the
            accepted shapes.
        method : str or callable, optional
            Quadrature rule name (``"simpson"``, ``"trapezoid"``,
            ``"gauss-legendre"``) or a callable returning
            ``(coords_1d, weights_1d)`` for ``(a, b, N)``.
        N : int, optional
            Number of quadrature points per dimension (odd for Simpson).
        chunk_size_params : int, optional
            Number of parameter points processed in one chunk.  Reducing this
            lowers peak memory usage at the cost of more Python loops.
        chunk_size_points : int or None, optional
            Maximum number of grid points processed per chunk.  ``None`` means
            "all at once".
        device : torch.device or None, optional
            Device used for all internal tensors.  If ``None``, CUDA is used
            when available, otherwise CPU.
        dtype : torch.dtype or None, optional
            Floating dtype for internal computations.  Defaults to ``float32``
            on CUDA and ``float64`` on CPU.  A complex dtype is accepted (the
            parameters are then complex and the quadrature nodes use the
            matching real dtype).

        Returns
        -------
        re, im : torch.Tensor
            Tensors with the same batch shape as the input parameters,
            containing the real and imaginary parts of the integral.

        USAGE
        =====
        1. 1D integral without parameters (single value)::

            ````python
            import torchsympy
            from sympy import symbols, Integral, exp, oo

            x = symbols("x", real=True)
            lt = torchsympy.TorchSymPy()
            texpr = lt.torchify(Integral(exp(-x**2), (x, -oo, oo)))

            re, im = texpr.torch_integrate_batched(params_values=None, N=121)
            ````

        2. 1D integral evaluated on a grid of parameter values::

            ````python
            import torchsympy, torch
            from sympy import symbols, Integral, exp, I, oo

            x, k = symbols("x k", real=True)
            lt = torchsympy.TorchSymPy()
            texpr = lt.torchify(Integral(exp(-x**2) * exp(I * k * x), (x, -oo, oo)))

            k_grid = torch.linspace(-5.0, 5.0, 250).unsqueeze(-1)  # (250, 1)
            re, im = texpr.torch_integrate_batched(
                params_values=k_grid, N=51,
                chunk_size_params=64, chunk_size_points=10_000,
            )
            ````
        """
        if self.dim <= 0:
            raise ValueError(f"self.dim must be positive; got dim={self.dim}")
        if self.n_params < 0:
            raise ValueError(f"self.n_params must be non-negative; got n_params={self.n_params}")
        if len(self.domain) != self.dim:
            raise ValueError(f"self.domain must have length dim={self.dim}; got {len(self.domain)}")
        if chunk_size_params <= 0:
            raise ValueError(f"chunk_size_params must be positive; got {chunk_size_params}")
        if chunk_size_points is not None and chunk_size_points <= 0:
            raise ValueError(f"chunk_size_points must be positive; got {chunk_size_points}")

        device, compute_dtype, param_dtype = _resolve_device_dtype(device, dtype)
        complex_dtype = _complex_dtype(compute_dtype)

        params_flat, batch_shape = self._normalize_params_values(
            params_values, self.n_params, device=device, dtype=param_dtype
        )
        batch_size = int(params_flat.shape[0])

        rule = self._resolve_quadrature_rule(method)
        rule_key = method if callable(method) else str(method).strip().lower()

        coords_1d: List[torch.Tensor] = []
        weights_1d: List[torch.Tensor] = []
        for bounds in self.domain:
            lower_bound, upper_bound = (float(bounds[0]), float(bounds[1]))
            coords, weights = rule(
                lower_bound, upper_bound, N, device=device, dtype=compute_dtype
            )
            coords = torch.as_tensor(coords, device=device, dtype=compute_dtype).reshape(-1)
            weights = torch.as_tensor(weights, device=device, dtype=compute_dtype).reshape(-1)
            if coords.numel() != weights.numel():
                raise ValueError(
                    "Quadrature rule returned mismatched coordinates/weights lengths: "
                    f"{coords.numel()} vs {weights.numel()}"
                )
            coords_1d.append(coords)
            weights_1d.append(weights)

        grid_key = (
            rule_key,
            tuple((float(bounds[0]), float(bounds[1])) for bounds in self.domain),
            int(N),
            device,
            compute_dtype,
        )
        n_points, take_grid = self._cached_grid_chunker(grid_key, coords_1d, weights_1d)
        points_per_chunk = n_points if chunk_size_points is None else chunk_size_points
        logger.debug(
            "batched quadrature: dim={} n_points={} batch={} dtype={} device={}",
            self.dim, n_points, batch_size, compute_dtype, device,
        )

        re_out = torch.empty(batch_size, device=device, dtype=compute_dtype)
        im_out = torch.zeros(batch_size, device=device, dtype=compute_dtype)

        for start_param in range(0, batch_size, chunk_size_params):
            stop_param = min(batch_size, start_param + chunk_size_params)
            n_batch = stop_param - start_param
            param_chunk = params_flat[start_param:stop_param, :]
            param_args = [param_chunk[:, index].unsqueeze(0) for index in range(self.n_params)]

            re_acc = torch.zeros(n_batch, device=device, dtype=compute_dtype)
            im_acc = torch.zeros(n_batch, device=device, dtype=compute_dtype)
            saw_complex = False

            for start_point in range(0, n_points, points_per_chunk):
                stop_point = min(n_points, start_point + points_per_chunk)
                grid_chunk, weight_chunk = take_grid(start_point, stop_point)
                n_pts = stop_point - start_point

                var_args = [grid_chunk[:, index].unsqueeze(1) for index in range(self.dim)]
                values = self._align_values(
                    self.func(*var_args, *param_args), n_pts, n_batch, device=device
                )

                # A single mat-vec is cheaper (time and memory) than an
                # elementwise product followed by a sum.
                if torch.is_complex(values):
                    saw_complex = True
                    contribution = weight_chunk.to(complex_dtype) @ values.to(complex_dtype)
                    re_acc += contribution.real
                    im_acc += contribution.imag
                else:
                    re_acc += weight_chunk @ values.to(compute_dtype)

            re_out[start_param:stop_param] = re_acc
            if saw_complex:
                im_out[start_param:stop_param] = im_acc

        return re_out.reshape(batch_shape), im_out.reshape(batch_shape)

    def torch_differentiate(self, domain_points=None, params_values=None, argnums=None):
        """Batched multi-dimensional differentiation of this ``TorchExpr``.

        Computes the batched Jacobian of the underlying PyTorch callable
        (``self.func``) over a multi-dimensional grid using ``torch.func``'s
        ``vmap`` and ``jacrev``.  It works for any function represented by this
        ``TorchExpr`` (an integrand, or a plain function wrapped through
        :meth:`TorchSymPy.torchify_callable`).

        Parameters
        ----------
        domain_points : iterable of torch.Tensor or None
            Variable tensors passed to the function (``self.dim`` of them).
        params_values : iterable of torch.Tensor or None
            Parameter tensors passed to the function (``self.n_params`` of them).
        argnums : int or tuple of ints, optional
            Indices of the arguments to differentiate with respect to.  If
            ``None``, differentiates with respect to all parameters (or all
            variables when ``n_params == 0``).  Indices ``0 .. dim-1`` are the
            domain points, ``dim .. dim+n_params-1`` the parameters.

        Returns
        -------
        torch.Tensor or tuple of torch.Tensor
            The batched Jacobian.  A tuple is returned exactly when ``argnums``
            is a tuple/list (including a tuple of length one), matching
            ``torch.func.jacrev``.
        """
        from torch.func import jacrev, vmap

        points = list(domain_points) if domain_points is not None else []
        params = list(params_values) if params_values is not None else []

        if len(points) != self.dim:
            raise ValueError(f"Expected {self.dim} domain_points, got {len(points)}")
        if len(params) != self.n_params:
            raise ValueError(f"Expected {self.n_params} params_values, got {len(params)}")

        all_args = points + params
        if not all_args:
            raise ValueError("No arguments to differentiate.")

        if argnums is None:
            targets = (
                tuple(range(self.dim, self.dim + self.n_params))
                if self.n_params > 0
                else tuple(range(self.dim))
            )
            # A single target yields a plain tensor rather than a 1-tuple.
            argnums = targets[0] if len(targets) == 1 else targets
        # ``jacrev`` returns a tuple whenever ``argnums`` is a sequence — even
        # for a sequence of length one.
        returns_tuple = isinstance(argnums, (list, tuple))
        if returns_tuple:
            argnums = tuple(argnums)

        # Differentiation needs floating point inputs; broadcast to a common shape.
        tensor_args = [_as_float_tensor(arg) for arg in all_args]
        broadcast_args = torch.broadcast_tensors(*tensor_args)
        batch_shape = broadcast_args[0].shape

        jacobian_fn = jacrev(self.func, argnums=argnums)

        if not batch_shape:  # scalar inputs
            return jacobian_fn(*broadcast_args)

        flat_args = [arg.reshape(-1) for arg in broadcast_args]
        flat_result = vmap(jacobian_fn)(*flat_args)

        if returns_tuple:
            return tuple(
                part.reshape(batch_shape + part.shape[1:]) for part in flat_result
            )
        return flat_result.reshape(batch_shape + flat_result.shape[1:])


class TorchSymPy:
    """Compiler from SymPy expressions/integrals to :class:`TorchExpr` objects."""

    #: memoises ``lambdify`` for repeated identical expressions built with the
    #: default module mapping.
    _lambdify_cache: Dict[Tuple[Any, Any, Tuple[Any, ...]], Callable] = {}
    _LAMBDIFY_CACHE_MAXSIZE = 256

    def __init__(
        self,
        *,
        separable_integration: bool = True,
        separable_tol: float = SEPARABLE_TOL,
        separable_start_nodes: int = SEPARABLE_START_NODES,
        separable_max_nodes: int = SEPARABLE_MAX_NODES,
        shift_invariant_integration: bool = True,
        shift_invariant_tol: float = SHIFT_INVARIANT_TOL,
        shift_invariant_start_nodes: int = SHIFT_INVARIANT_START_NODES,
        shift_invariant_max_nodes: int = SHIFT_INVARIANT_MAX_NODES,
    ) -> None:
        """Create a compiler.

        Parameters
        ----------
        separable_integration : bool, optional
            When ``True`` (the default) :meth:`eval_numeric` first tries to
            factorise every multi-dimensional ``Integral`` into independent 1D
            integrals and to refine each of them until convergence.  This is
            what makes oscillatory (diffraction) integrands correct at short
            propagation distances; see the module-level comment above
            :func:`_separable_plan`.  Set it to ``False`` to restore the
            previous behaviour exactly (fixed-``N`` tensor-product quadrature).
        separable_tol : float, optional
            Relative tolerance of the refinement loop.
        separable_start_nodes, separable_max_nodes : int, optional
            First and last node count per axis of the refinement loop.
        shift_invariant_integration : bool, optional
            When ``True`` (the default) an ``Integral`` that is *not* separable
            but only depends on the differences between its integration
            variables and the parameters - the Rayleigh-Sommerfeld kernel is
            the typical case - is evaluated through one cumulative integral of
            the kernel that is shared by every parameter value, again refined
            until it converges; see the comment above
            :func:`_shift_invariant_plan`.  Set it to ``False`` to restore the
            previous behaviour exactly.
        shift_invariant_tol : float, optional
            Relative tolerance of that refinement loop.
        shift_invariant_start_nodes, shift_invariant_max_nodes : int, optional
            First and last Gauss-Legendre node count per panel.
        """
        self.separable_integration = bool(separable_integration)
        self.separable_tol = float(separable_tol)
        self.separable_start_nodes = int(separable_start_nodes)
        self.separable_max_nodes = int(separable_max_nodes)
        self.shift_invariant_integration = bool(shift_invariant_integration)
        self.shift_invariant_tol = float(shift_invariant_tol)
        self.shift_invariant_start_nodes = int(shift_invariant_start_nodes)
        self.shift_invariant_max_nodes = int(shift_invariant_max_nodes)

    # ------------------------------------------------------------------
    # Change of variables
    # ------------------------------------------------------------------
    _COV_ALIASES = {
        "tangent": "tangent",
        "tan": "tangent",
        "algebraic": "algebraic",
        "rational": "algebraic",
        "tanh-sinh": "tanh-sinh",
        "tanh_sinh": "tanh-sinh",
        "tanhsinh": "tanh-sinh",
    }

    @staticmethod
    def _normalize_cov_method(change_of_variables_method: str) -> str:
        """Canonical name of a change-of-variables method."""
        method = str(change_of_variables_method).strip().lower()
        try:
            return TorchSymPy._COV_ALIASES[method]
        except KeyError:
            raise ValueError(
                "Unknown change_of_variables_method. Expected one of "
                "{'tangent', 'algebraic', 'tanh-sinh'}. "
                f"Got: {change_of_variables_method}"
            ) from None

    @staticmethod
    def _inset_open_interval(lower: float, upper: float, eps: float) -> List[float]:
        """Shrink an open interval by ``eps`` on both sides.

        ``eps`` is clipped so that the resulting interval never degenerates.
        """
        lower, upper = float(lower), float(upper)
        eps = min(abs(float(eps)), 0.25 * (upper - lower))
        return [lower + eps, upper - eps]

    @staticmethod
    def _is_finite_number(bound) -> bool:
        """``True`` when ``bound`` is a concrete finite number."""
        if bound in (oo, -oo):
            return False
        try:
            return math.isfinite(float(bound))
        except (TypeError, ValueError):
            return False

    def _transform_limit_with_method(
        self, variable, lower, upper, *, change_of_variables_method: str, eps: float
    ):
        """Return ``(new_variable, substitution, jacobian, finite_domain)``.

        Three cases are handled:

        * numeric finite limits — no transformation;
        * *symbolic* finite limits (e.g. a nested integral ``(y, 0, x)``) — an
          affine map onto ``[0, 1]``, which keeps the symbolic bound as a
          parameter;
        * infinite / semi-infinite limits — the requested change of variables.
        """
        method = self._normalize_cov_method(change_of_variables_method)

        lower_is_inf = lower == -oo or lower == oo
        upper_is_inf = upper == oo or upper == -oo

        if not lower_is_inf and not upper_is_inf:
            if self._is_finite_number(lower) and self._is_finite_number(upper):
                return variable, variable, Integer(1), [float(lower), float(upper)]

            # Symbolic finite bounds: u in [0, 1] with x = lower + (upper-lower)*u.
            u_var = Symbol(f"u_{variable.name}", real=True)
            width = upper - lower
            mapped = lower + width * u_var
            logger.debug("affine change of variables for symbolic limits of {}", variable)
            return u_var, mapped, width, [0.0, 1.0]

        t_var = Symbol(f"t_{variable.name}", real=True)

        if method == "tangent":
            if lower == -oo and upper == oo:
                mapped = tan(t_var)
                open_interval = (float(-pi / 2), float(pi / 2))
            elif upper == oo:
                mapped = lower + tan(t_var) ** 2
                open_interval = (0.0, float(pi / 2))
            elif lower == -oo:
                mapped = upper - tan(t_var) ** 2
                open_interval = (0.0, float(pi / 2))
            else:
                raise ValueError(f"Unexpected limit pattern in tangent transform: ({lower}, {upper})")

        elif method == "algebraic":
            if lower == -oo and upper == oo:
                mapped = t_var / (1 - t_var ** 2)
                open_interval = (-1.0, 1.0)
            elif upper == oo:
                mapped = lower + t_var / (1 - t_var)
                open_interval = (0.0, 1.0)
            elif lower == -oo:
                mapped = upper - t_var / (1 - t_var)
                open_interval = (0.0, 1.0)
            else:
                raise ValueError(f"Unexpected limit pattern in algebraic transform: ({lower}, {upper})")

        else:  # "tanh-sinh"
            # Finite-interval parameter t in (-1,1) or (0,1) mapped through sinh(atanh(t)).
            core = sinh(atanh(t_var))
            if lower == -oo and upper == oo:
                mapped = core
                open_interval = (-1.0, 1.0)
            elif upper == oo:
                mapped = lower + core ** 2
                open_interval = (0.0, 1.0)
            elif lower == -oo:
                mapped = upper - core ** 2
                open_interval = (0.0, 1.0)
            else:
                raise ValueError(f"Unexpected limit pattern in tanh-sinh transform: ({lower}, {upper})")

        jacobian = diff(mapped, t_var)
        domain = self._inset_open_interval(open_interval[0], open_interval[1], eps)
        return t_var, mapped, jacobian, domain

    # ------------------------------------------------------------------
    # lambdify helpers
    # ------------------------------------------------------------------
    def _default_modules(self) -> List[Dict[str, Any]]:
        """Default SymPy → torch mapping for ``lambdify``.

        The returned list is passed directly to SymPy's ``lambdify`` as the
        ``modules`` argument.  It exposes a subset of functions implemented
        with torch operations so that the resulting numerical function is
        differentiable and GPU-friendly.

        Subclass ``TorchSymPy`` and override this method to customise it.
        """
        return [_default_module_mapping()]

    def _merge_lambdify_modules(self, modules, extra_mapping):
        """Merge ``modules`` with ``extra_mapping`` at the highest priority."""
        if not extra_mapping:
            return modules
        if modules is None:
            modules = self._default_modules()
        if isinstance(modules, (list, tuple)):
            return [dict(extra_mapping), *modules]
        return [dict(extra_mapping), modules]

    def _lambdify(self, args, expr, modules, *, cacheable: bool):
        """``lambdify`` with an optional memo for repeated identical calls."""
        if not cacheable:
            return lambdify(args, expr, modules=modules)

        # SymPy objects hash structurally (and cache their hash), which makes
        # this much cheaper than serialising the expression.  The class is part
        # of the key so that a subclass overriding ``_default_modules`` never
        # shares compiled functions with the base class.
        key = (type(self), expr, tuple(args))
        cached = self._lambdify_cache.get(key)
        if cached is None:
            cached = lambdify(args, expr, modules=modules)
            if len(self._lambdify_cache) >= self._LAMBDIFY_CACHE_MAXSIZE:
                self._lambdify_cache.clear()
            self._lambdify_cache[key] = cached
        else:
            logger.debug("lambdify cache hit")
        return cached

    # ------------------------------------------------------------------
    # Nested definite integrals
    # ------------------------------------------------------------------
    def _build_nested_definite_integral_callable(
        self, texpr: TorchExpr, n_params: int, *, inner_N: int = 61
    ):
        """Numeric callable for ``lambdify`` evaluating a nested definite integral.

        The whole batch of outer evaluation points is integrated in one batched
        call instead of one torchquad call per point.
        """

        def eval_nested(*args):
            if len(args) != n_params:
                raise ValueError(f"Expected {n_params} nested integral params; got {len(args)}")

            def combine(re_val, im_val):
                # Stay real when the inner integral is: a complex outer
                # integrand costs twice as much to evaluate.
                if bool(torch.any(im_val != 0)):
                    return re_val + 1j * im_val
                return re_val

            if n_params == 0:
                return combine(*texpr.torch_integrate_batched(N=inner_N))

            tensors = torch.broadcast_tensors(*[_as_float_tensor(arg) for arg in args])
            batch_shape = tuple(tensors[0].shape)
            device = tensors[0].device
            dtype = torch.promote_types(
                tensors[0].dtype, torch.float32 if device.type == "cuda" else torch.float64
            )

            params = torch.stack([t.reshape(-1) for t in tensors], dim=-1)
            re_val, im_val = texpr.torch_integrate_batched(
                params_values=params, N=inner_N, device=device, dtype=dtype
            )
            return combine(re_val, im_val).reshape(batch_shape)

        return eval_nested

    def _replace_nested_definite_integrals(
        self,
        expr,
        *,
        inner_N: int = 61,
        change_of_variables_method: str = "tangent",
        cov_eps: float = 1e-7,
    ):
        """Replace nested definite ``Integral`` nodes with numeric callables.

        Any nested integral with a non-definite limit triggers a ``ValueError``.
        Symbolic (outer-variable dependent) bounds are supported.
        """
        if not hasattr(expr, "has") or not expr.has(Integral):
            return expr, {}

        expr_work = expr
        nested_mapping: Dict[str, Callable] = {}
        counter = 0

        while expr_work.has(Integral):
            leaf_integrals = [
                integral
                for integral in expr_work.atoms(Integral)
                if not integral.function.has(Integral)
            ]
            if not leaf_integrals:
                break

            for inner_integral in leaf_integrals:
                for limit in inner_integral.limits:
                    if len(limit) != 3:
                        raise ValueError(
                            "The integrand must not contain unevaluated indefinite integrals. "
                            f"Found nested integral with non-definite limit: {inner_integral}"
                        )

                integration_vars = {limit[0] for limit in inner_integral.limits}
                param_symbols = sorted(
                    inner_integral.free_symbols - integration_vars, key=lambda symbol: symbol.name
                )

                nested_texpr = self.torchify(
                    inner_integral,
                    params=param_symbols,
                    change_of_variables_method=change_of_variables_method,
                    cov_eps=cov_eps,
                )
                nested_name = f"_nested_definite_integral_{counter}"
                counter += 1
                logger.debug("replacing nested integral {} by {}", inner_integral, nested_name)

                nested_mapping[nested_name] = self._build_nested_definite_integral_callable(
                    nested_texpr, len(param_symbols), inner_N=inner_N
                )
                expr_work = expr_work.xreplace(
                    {inner_integral: Function(nested_name)(*param_symbols)}
                )

        if expr_work.has(Integral):
            raise ValueError(
                "The integrand must not contain unevaluated integrals. "
                "Only nested definite integrals are supported."
            )

        return expr_work, nested_mapping

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    @staticmethod
    def _absorb_scalar_mul_into_integral(expr_candidate):
        """Rewrite ``prefactor * Integral(...)`` as ``Integral(prefactor*integrand, ...)``.

        This is a convenience for common patterns like ``Integral(f, ...)/(2*pi)``
        which SymPy stores as ``Mul(Integral(...), 1/(2*pi))``.  The rewrite is
        conservative: the prefactor is only absorbed when it is independent of
        the integral's dummy variables.
        """
        if not isinstance(expr_candidate, Mul):
            return expr_candidate

        integral_factors = [arg for arg in expr_candidate.args if isinstance(arg, Integral)]
        if len(integral_factors) != 1:
            return expr_candidate

        integral_factor = integral_factors[0]
        prefactor = Mul(*[arg for arg in expr_candidate.args if arg is not integral_factor])

        dummy_vars = {limit[0] for limit in integral_factor.limits}
        if not prefactor.free_symbols.isdisjoint(dummy_vars):
            return expr_candidate

        return Integral(prefactor * integral_factor.function, *integral_factor.limits)

    @classmethod
    def _extract_integral(cls, expr):
        """Return the ``Integral`` carried by ``expr``, or ``None``."""
        if isinstance(expr, Integral):
            return expr
        if isinstance(expr, Eq):
            if isinstance(expr.rhs, Integral):
                return expr.rhs
            candidate = cls._absorb_scalar_mul_into_integral(expr.rhs)
            return candidate if isinstance(candidate, Integral) else None
        candidate = cls._absorb_scalar_mul_into_integral(expr)
        return candidate if isinstance(candidate, Integral) else None

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
        """Convert a SymPy object into a torch-compatible callable or ``TorchExpr``.

        The change of variables for infinite / semi-infinite limits is done at
        the SymPy level (once), then the transformed expression is
        ``lambdify``'d.

        ``expr`` can be:

        - a plain SymPy expression (e.g. ``exp(-x**2)``),
        - a SymPy ``Integral`` (``Integral(f(x), (x, -oo, oo))``), or
        - a SymPy ``Eq`` whose right-hand side is an ``Integral``.

        If an ``Integral`` (or an equation with an integral) is passed, the
        integrand and limits are extracted internally, so there is no need to
        use ``.function`` / ``.limits`` manually.

        Parameters
        ----------
        expr : sympy.Expr | sympy.Integral | sympy.Eq
            Symbolic expression (may be complex-valued) or an integral.
        variables : list[sympy.Symbol] or None
            Integration variables (order matters).  If ``expr`` is an
            ``Integral`` and ``variables`` is ``None``, they are inferred from
            the integral limits.
        limits : list[tuple] or None
            ``(var, lower, upper)`` for each variable.  ``oo`` / ``-oo``
            trigger an automatic change of variables; finite symbolic bounds
            trigger an affine one.  If ``expr`` is an ``Integral`` and
            ``limits`` is ``None``, they are taken from ``expr.limits``.
        params : list[sympy.Symbol] or None
            Extra symbolic parameters that appear in the transformed
            expression but are **not** integrated over.  If ``None``, they are
            inferred from the free symbols (sorted by name).
        modules : list[dict] or None
            Custom ``lambdify`` module list.  ``None`` → default torch mapping
            from :meth:`_default_modules`.
        change_of_variables_method : str, optional
            ``"tangent"``, ``"algebraic"`` or ``"tanh-sinh"``.
        cov_eps : float, optional
            Small inward shift for transformed open intervals to avoid
            evaluating at a singular endpoint.

        Returns
        -------
        TorchExpr
            If limits are provided (directly or via an ``Integral``).
        callable
            If ``limits`` is ``None``: a plain torch-ready function.

        USAGE
        =====
        1. Plain function (no integration)::

            ````python
            import torchsympy
            from sympy import symbols, exp

            x = symbols("x", real=True)
            lt = torchsympy.TorchSymPy()
            f = lt.torchify(exp(-x**2), variables=[x])
            ````

        2. Integral with finite limits::

            ````python
            texpr = lt.torchify(Integral(exp(-x**2), (x, 0, 1)))
            ````

        3. Integral with infinite limits (automatic change of variables)::

            ````python
            texpr = lt.torchify(Integral(exp(-x**2), (x, -oo, oo)))
            ````

        4. Equation with an integral on the right-hand side::

            ````python
            texpr = lt.torchify(Eq(y, Integral(exp(-x**2), (x, -oo, oo))))
            ````
        """
        integral = self._extract_integral(expr)
        nested_modules: Dict[str, Callable] = {}

        if integral is not None:
            base_expr = integral.function
            if base_expr.has(Integral):
                # Evaluate nested definite integrals numerically.
                base_expr, nested_modules = self._replace_nested_definite_integrals(
                    base_expr,
                    change_of_variables_method=change_of_variables_method,
                    cov_eps=cov_eps,
                )
            if limits is None:
                limits = list(integral.limits)
            if variables is None:
                variables = [limit[0] for limit in limits]
        else:
            base_expr = expr

        if variables is None:
            raise ValueError("'variables' must be provided when expr is not an Integral")

        use_default_modules = modules is None
        if use_default_modules:
            modules = self._default_modules()
        modules = self._merge_lambdify_modules(modules, nested_modules)
        # Only memoise compilations that use the (stateless) default mapping.
        cacheable = use_default_modules and not nested_modules

        variables = list(variables)

        # ------------------------------------------------------------------
        # Simple case: no limits → plain lambdify
        # ------------------------------------------------------------------
        if limits is None:
            return self._lambdify(tuple(variables), base_expr, modules, cacheable=cacheable)

        # ------------------------------------------------------------------
        # With limits: symbolic change of variables (done once, fast forever)
        # ------------------------------------------------------------------
        for limit in limits:
            if len(limit) != 3:
                raise ValueError(f"Expected definite limits (var, lower, upper); got {limit}")
        limit_map = {limit[0]: (limit[1], limit[2]) for limit in limits}
        missing = [v for v in variables if v not in limit_map]
        if missing:
            raise ValueError(f"No integration limits given for {missing}")

        new_vars: List[Symbol] = []
        domain: List[List[float]] = []
        substitutions: List[Tuple[Symbol, Any]] = []
        jacobian = Integer(1)

        for variable in variables:
            lower, upper = limit_map[variable]
            new_var, mapped_expr, jacobian_expr, domain_entry = self._transform_limit_with_method(
                variable,
                lower,
                upper,
                change_of_variables_method=change_of_variables_method,
                eps=cov_eps,
            )
            new_vars.append(new_var)
            domain.append(domain_entry)
            if new_var != variable:
                substitutions.append((variable, mapped_expr))
            jacobian = jacobian * jacobian_expr

        # Substitute in the order of ``variables`` (innermost limit first), so
        # that a limit depending on an outer variable — e.g. ``(y, 0, x)`` —
        # gets that outer variable transformed as well.  Passing the whole dict
        # to ``subs`` would leave the order up to SymPy.
        expr_work = base_expr * jacobian
        for variable, mapped_expr in substitutions:
            expr_work = expr_work.subs(variable, mapped_expr)

        if params is None:
            params = sorted(expr_work.free_symbols - set(new_vars), key=lambda s: s.name)
        else:
            params = list(params)

        # Lambdify: new_vars first, then params.
        arglist = tuple(new_vars) + tuple(params)
        func = self._lambdify(arglist, expr_work, modules, cacheable=cacheable)
        logger.debug("torchified {} with dim={} n_params={}", integral or expr, len(new_vars), len(params))

        return TorchExpr(
            func=func,
            domain=domain,
            dim=len(new_vars),
            n_params=len(params),
            sympy_integrand=expr_work,
            sympy_integral=integral,
            variables=new_vars,
        )

    def torchify_callable(
        self, integrand: Callable, domain: List[List[float]], n_params: int = 0
    ) -> TorchExpr:
        """Create a ``TorchExpr`` directly from a Python callable (no SymPy).

        This hooks a custom PyTorch function into the ``TorchExpr``
        integration backends without any symbolic processing.

        The callable must accept the integration variables first (as many as
        ``domain`` has entries), followed by ``n_params`` parameters.

        Parameters
        ----------
        integrand : Callable
            ``f(var_0, ..., var_{dim-1}, param_0, ..., param_{n_params-1})``.
        domain : list[list[float]]
            Finite integration box, e.g. ``[[0.0, 1.0], [-5.0, 5.0]]``.  Its
            length sets the integration dimensionality.
        n_params : int, optional
            Number of extra parameters accepted after the variables.

        Returns
        -------
        TorchExpr
            Ready for fast numerical integration.

        USAGE
        =====
        1. Simple 2D integration with parameters::

            ````python
            import torch, torchsympy
            lt = torchsympy.TorchSymPy()

            def my_integrand(x, y, alpha, beta):
                return torch.exp(-alpha * x**2 - beta * y**2)

            texpr = lt.torchify_callable(
                integrand=my_integrand,
                domain=[[-5.0, 5.0], [-5.0, 5.0]],  # 2 variables: x, y
                n_params=2,                          # 2 parameters: alpha, beta
            )
            re, im = texpr.torchquad_integrate(params_values=[1.0, 0.5])
            ````
        """
        if not callable(integrand):
            raise TypeError("integrand must be callable")
        if n_params < 0:
            raise ValueError(f"n_params must be non-negative; got {n_params}")
        domain = [[float(a), float(b)] for a, b in domain]
        if not domain:
            raise ValueError("domain must contain at least one [a, b] pair")

        return TorchExpr(
            func=integrand,
            domain=domain,
            dim=len(domain),
            n_params=n_params,
            sympy_integrand=None,
            sympy_integral=None,
            variables=None,
        )

    # ------------------------------------------------------------------
    # Composite expressions containing Integral nodes
    # ------------------------------------------------------------------
    @staticmethod
    def _split_params_for_lambdify(params_values, n_params: int) -> List[Any]:
        """Return one entry per symbolic parameter.

        Accepts the same shapes as the solvers: ``None``, a sequence with one
        entry per parameter, or a tensor/array of shape ``(..., n_params)``.
        """
        if params_values is None:
            return []
        if torch.is_tensor(params_values) or isinstance(params_values, np.ndarray):
            tensor = torch.as_tensor(params_values)
            if tensor.ndim > 0 and tensor.shape[-1] == n_params:
                return [tensor[..., index] for index in range(n_params)]
            if n_params == 1:
                return [tensor]
            raise ValueError(
                f"params_values must have last dimension n_params={n_params}; "
                f"got shape {tuple(tensor.shape)}"
            )
        if not isinstance(params_values, (list, tuple)):
            params_values = [params_values]
        params_list = list(params_values)
        if len(params_list) != n_params:
            raise ValueError(f"Expected {n_params} params; got {len(params_list)}")
        return params_list

    # ------------------------------------------------------------------
    # Separable evaluation of oscillatory integrals
    # ------------------------------------------------------------------
    @staticmethod
    def _quadrature_1d(func, lower, upper, param_columns, n_nodes, *, device, dtype):
        """Composite Simpson rule for ``func(variable, *params)`` on one axis.

        ``param_columns`` holds one 1-D tensor per parameter of ``func``; the
        rule is applied to all of their entries at once, in chunks of at most
        :data:`MAX_ELEMENTS_PER_EVALUATION` integrand values.
        """
        coords, weights = _cached_rule("simpson", lower, upper, n_nodes, device, dtype)
        n_values = int(param_columns[0].numel()) if param_columns else 1
        budget = max(1, MAX_ELEMENTS_PER_EVALUATION // n_nodes)

        pieces: List[torch.Tensor] = []
        for start in range(0, n_values, budget):
            stop = min(n_values, start + budget)
            args = [column[start:stop].unsqueeze(0) for column in param_columns]
            values = func(coords.unsqueeze(1), *args)
            if not torch.is_tensor(values):
                values = torch.as_tensor(values, device=device, dtype=dtype)
            values = torch.broadcast_to(values, (n_nodes, stop - start))
            pieces.append((weights.unsqueeze(1) * values).sum(dim=0))
        return torch.cat(pieces)

    def _quadrature_1d_refined(self, func, lower, upper, param_columns, *, device, dtype, label):
        """Apply :meth:`_quadrature_1d` with doubling until the result settles.

        The number of oscillations of a diffraction integrand depends on the
        wavelength, the propagation distance and the screen coordinate, so a
        fixed node count is either wasteful or aliased.  Doubling until two
        successive rules agree to ``separable_tol`` picks the resolution the
        integrand actually needs.
        """
        n_nodes = max(3, int(self.separable_start_nodes))
        if n_nodes % 2 == 0:
            n_nodes += 1
        previous: Optional[torch.Tensor] = None

        while True:
            value = self._quadrature_1d(
                func, lower, upper, param_columns, n_nodes, device=device, dtype=dtype
            )
            if previous is not None:
                scale = torch.max(torch.abs(value))
                difference = torch.max(torch.abs(value - previous))
                if bool(difference <= self.separable_tol * torch.clamp(scale, min=1e-300)):
                    logger.debug(
                        "separable quadrature over {} converged with {} nodes", label, n_nodes
                    )
                    return value
            if n_nodes >= self.separable_max_nodes:
                logger.warning(
                    "separable quadrature over {} still changing at {} nodes per axis "
                    "(separable_max_nodes); the result may be under-resolved",
                    label, n_nodes,
                )
                return value
            previous = value
            n_nodes = 2 * (n_nodes - 1) + 1

    def _evaluate_separable_integral(
        self,
        integral,
        free_symbols,
        params_list,
        *,
        dtype=None,
        device=None,
        **unused_solver_kwargs,
    ):
        """Evaluate ``integral`` as a product of refined 1D integrals.

        Returns ``None`` when the integrand is not separable (or the parameters
        are not plain real tensors), so that the caller can fall back to the
        tensor-product solvers.
        """
        if not getattr(self, "separable_integration", True):
            return None
        plan = _separable_plan(integral)
        if plan is None:
            return None

        device, compute_dtype, _ = _resolve_device_dtype(device, dtype)

        tensors: Dict[Symbol, torch.Tensor] = {}
        for symbol, value in zip(free_symbols, params_list):
            tensor = value if torch.is_tensor(value) else torch.as_tensor(np.asarray(value))
            if torch.is_complex(tensor):
                return None
            tensors[symbol] = tensor.to(device=device, dtype=compute_dtype)

        shapes = [tuple(tensor.shape) for tensor in tensors.values() if tensor.ndim > 0]
        batch_shape = tuple(torch.broadcast_shapes(*shapes)) if shapes else ()
        n_values = int(np.prod(batch_shape)) if batch_shape else 1
        flat = {
            symbol: torch.broadcast_to(tensor, batch_shape).reshape(-1)
            if batch_shape
            else tensor.reshape(1)
            for symbol, tensor in tensors.items()
        }

        def lambdified(expression, arguments):
            return self._lambdify(
                tuple(arguments), expression, self._default_modules(), cacheable=True
            )

        total = torch.ones(1, device=device, dtype=compute_dtype)
        for factor in plan.factors:
            used = sorted(
                (factor.amplitude.free_symbols | factor.phase.free_symbols) - {factor.variable},
                key=lambda symbol: symbol.name,
            )
            if any(symbol not in flat for symbol in used):
                return None

            integrand = factor.amplitude
            if factor.phase != 0:
                integrand = integrand * sympy_exp(sympy_I * factor.phase)
            func = lambdified(integrand, (factor.variable, *used))

            if used:
                stacked = torch.stack([flat[symbol] for symbol in used], dim=1)
                unique, inverse = torch.unique(stacked, dim=0, return_inverse=True)
                columns = [unique[:, index].contiguous() for index in range(len(used))]
            else:
                columns, inverse = [], None

            value = self._quadrature_1d_refined(
                func,
                factor.lower,
                factor.upper,
                columns,
                device=device,
                dtype=compute_dtype,
                label=factor.variable,
            )
            if inverse is not None:
                value = value[inverse]
            total = total * value

        coefficient = plan.coeff
        if plan.const_phase != 0:
            coefficient = coefficient * sympy_exp(sympy_I * plan.const_phase)
        if coefficient != 1:
            used = sorted(coefficient.free_symbols, key=lambda symbol: symbol.name)
            if any(symbol not in flat for symbol in used):
                return None
            value = lambdified(coefficient, tuple(used))(*[flat[symbol] for symbol in used])
            if not torch.is_tensor(value):
                value = torch.as_tensor(value, device=device)
            total = total * value

        if plan.trig == "cos":
            total = total.real if torch.is_complex(total) else total
        elif plan.trig == "sin":
            total = total.imag if torch.is_complex(total) else torch.zeros_like(total)
        if plan.outer == "re":
            total = total.real if torch.is_complex(total) else total
        elif plan.outer == "im":
            total = total.imag if torch.is_complex(total) else torch.zeros_like(total)

        logger.debug(
            "evaluated {} as {} separable 1D integrals", integral, len(plan.factors)
        )
        return torch.broadcast_to(total.reshape(-1), (n_values,)).reshape(batch_shape)

    # ------------------------------------------------------------------
    # Shift-invariant (convolution) evaluation of oscillatory integrals
    # ------------------------------------------------------------------
    @staticmethod
    def _shift_invariant_pass(func, boundaries, indices, scalars, *, n_nodes, device, dtype):
        """One pass of the cumulative-panel rule with ``n_nodes`` nodes per panel.

        ``boundaries`` holds the panel boundaries of each axis, ``indices`` the
        index of the lower/upper corner of every screen point in them.  The
        kernel is evaluated once on the whole panel grid; the panel integrals
        are accumulated into the cumulative integral ``Phi`` and every screen
        point is then read off as a difference of ``2**d`` of its values.
        """
        nodes, weights = [], []
        for boundary in boundaries:
            axis_nodes, axis_weights = _panel_nodes_weights(boundary, n_nodes)
            nodes.append(axis_nodes)
            weights.append(axis_weights)

        if len(boundaries) == 1:
            n_panels = nodes[0].shape[0]
            budget = max(1, MAX_ELEMENTS_PER_EVALUATION // n_nodes)
            pieces = []
            for start in range(0, n_panels, budget):
                stop = min(n_panels, start + budget)
                values = func(nodes[0][start:stop].reshape(-1), *scalars)
                values = _broadcast_values(
                    values, ((stop - start) * n_nodes,), device=device, dtype=dtype
                ).reshape(stop - start, n_nodes)
                pieces.append(
                    (weights[0][start:stop].to(values.dtype) * values).sum(dim=1)
                )
            panels_1d = torch.cat(pieces)
            cumulative = torch.cat(
                [
                    torch.zeros(1, device=device, dtype=panels_1d.dtype),
                    panels_1d.cumsum(dim=0),
                ]
            )
            lower, upper = indices[0]
            return cumulative[upper] - cumulative[lower]

        n_panels_u, n_panels_v = nodes[0].shape[0], nodes[1].shape[0]
        n_v = n_panels_v * n_nodes
        v_nodes = nodes[1].reshape(-1)
        budget = max(1, MAX_ELEMENTS_PER_EVALUATION // (n_v * n_nodes))
        strips = []
        for start in range(0, n_panels_u, budget):
            stop = min(n_panels_u, start + budget)
            u_nodes = nodes[0][start:stop].reshape(-1)
            values = func(u_nodes.unsqueeze(1), v_nodes.unsqueeze(0), *scalars)
            values = _broadcast_values(
                values, ((stop - start) * n_nodes, n_v), device=device, dtype=dtype
            ).reshape(stop - start, n_nodes, n_v)
            strips.append(
                torch.einsum(
                    "sp,spv->sv", weights[0][start:stop].to(values.dtype), values
                )
            )
        strip = torch.cat(strips)
        panels_2d = torch.einsum(
            "tq,stq->st",
            weights[1].to(strip.dtype),
            strip.reshape(n_panels_u, n_panels_v, n_nodes),
        )
        cumulative = torch.zeros(
            (n_panels_u + 1, n_panels_v + 1), device=device, dtype=panels_2d.dtype
        )
        cumulative[1:, 1:] = panels_2d.cumsum(dim=0).cumsum(dim=1)
        (lower_u, upper_u), (lower_v, upper_v) = indices
        return (
            cumulative[upper_u, upper_v]
            - cumulative[lower_u, upper_v]
            - cumulative[upper_u, lower_v]
            + cumulative[lower_u, lower_v]
        )

    def _evaluate_shift_invariant_integral(
        self,
        integral,
        free_symbols,
        params_list,
        *,
        dtype=None,
        device=None,
        **unused_solver_kwargs,
    ):
        """Evaluate a convolution integral through one cumulative panel rule.

        This is the path that keeps the Rayleigh-Sommerfeld integrand correct:
        its square root couples the two integration variables, so it is not
        separable, but it only depends on ``x - x0`` and ``y - y0``.  Returns
        ``None`` when the integrand is not of that form (or the parameters are
        not plain real tensors), so that the caller can fall back to the
        tensor-product solvers.
        """
        if not getattr(self, "shift_invariant_integration", True):
            return None
        plan = _shift_invariant_plan(integral, free_symbols)
        if plan is None:
            return None

        device, compute_dtype, _ = _resolve_device_dtype(device, dtype)

        tensors: Dict[Symbol, torch.Tensor] = {}
        for symbol, value in zip(free_symbols, params_list):
            tensor = value if torch.is_tensor(value) else torch.as_tensor(np.asarray(value))
            if torch.is_complex(tensor):
                return None
            tensors[symbol] = tensor.to(device=device, dtype=compute_dtype)
        if any(axis.parameter not in tensors for axis in plan.axes):
            return None

        shapes = [tuple(tensor.shape) for tensor in tensors.values() if tensor.ndim > 0]
        batch_shape = tuple(torch.broadcast_shapes(*shapes)) if shapes else ()
        n_values = int(np.prod(batch_shape)) if batch_shape else 1
        flat = {
            symbol: torch.broadcast_to(tensor, batch_shape).reshape(-1)
            if batch_shape
            else tensor.reshape(1)
            for symbol, tensor in tensors.items()
        }

        # Everything the kernel still depends on, apart from the difference
        # coordinates, has to be constant over the mesh - otherwise the single
        # cumulative integral below would not be shared by all screen points.
        coordinates = [axis.coordinate for axis in plan.axes]
        extra = sorted(
            plan.kernel.free_symbols - set(coordinates), key=lambda symbol: symbol.name
        )
        scalars = []
        for symbol in extra:
            if symbol not in flat:
                return None
            unique = torch.unique(flat[symbol])
            if unique.numel() != 1:
                return None
            scalars.append(unique.reshape(()))

        boundaries, indices, n_panels = [], [], 1
        for axis in plan.axes:
            parameter = flat[axis.parameter]
            corners = torch.cat(
                [
                    axis.sign * parameter + axis.offset_lower,
                    axis.sign * parameter + axis.offset_upper,
                ]
            )
            axis_boundaries, inverse = _panel_boundaries(
                corners.detach().cpu().numpy().astype(np.float64)
            )
            if axis_boundaries.size < 2:
                return None
            n_panels *= axis_boundaries.size - 1
            boundaries.append(
                torch.as_tensor(axis_boundaries, device=device, dtype=compute_dtype)
            )
            inverse_tensor = torch.as_tensor(
                np.asarray(inverse).reshape(-1), device=device, dtype=torch.long
            )
            indices.append((inverse_tensor[:n_values], inverse_tensor[n_values:]))

        func = self._lambdify(
            tuple(coordinates + extra), plan.kernel, self._default_modules(), cacheable=True
        )

        n_nodes = max(2, int(self.shift_invariant_start_nodes))
        previous: Optional[torch.Tensor] = None
        while True:
            value = self._shift_invariant_pass(
                func,
                boundaries,
                indices,
                scalars,
                n_nodes=n_nodes,
                device=device,
                dtype=compute_dtype,
            )
            if previous is not None:
                scale = torch.max(torch.abs(value))
                difference = torch.max(torch.abs(value - previous))
                if bool(
                    difference <= self.shift_invariant_tol * torch.clamp(scale, min=1e-300)
                ):
                    logger.debug(
                        "shift-invariant quadrature converged with {} nodes on each of "
                        "{} panels", n_nodes, n_panels,
                    )
                    break
            if n_nodes >= self.shift_invariant_max_nodes:
                logger.warning(
                    "shift-invariant quadrature still changing at {} nodes per panel "
                    "(shift_invariant_max_nodes); the result may be under-resolved",
                    n_nodes,
                )
                break
            previous = value
            n_nodes *= 2

        logger.debug("evaluated {} as a cumulative convolution integral", integral)
        return value.reshape(batch_shape) if batch_shape else value.reshape(())

    def eval_numeric(
        self,
        expr,
        params_values=None,
        solver: str = "batched",
        change_of_variables_method: str = "tangent",
        cov_eps: float = 1e-7,
        **solver_kwargs,
    ):
        """Evaluate a composite SymPy expression containing ``Integral`` terms.

        Walks the expression tree, finds all ``Integral`` nodes, torchifies and
        numerically integrates each one with the chosen solver, then assembles
        the final result.  For scalar parameters this returns a SymPy ``Float``;
        for batched tensor parameters it returns a ``torch.Tensor``.

        Parameters
        ----------
        expr : sympy.Expr | sympy.Eq
            Any SymPy expression (or ``Eq``) that may contain ``Integral``
            nodes at arbitrary depth, e.g.::

                C * Integral(sin(...), (x0, -1, 1), (y0, -1, 1))**2
                + C * Integral(cos(...), (x0, -1, 1), (y0, -1, 1))**2

        params_values : list | tuple | torch.Tensor | None, optional
            Values for the free symbols (parameters) of the integrals, ordered
            by symbol name:

            - ``None`` — no external parameters,
            - a list of scalars — single-point evaluation, returns a SymPy
              ``Float``,
            - a list of tensors / arrays, or a tensor of shape
              ``(..., n_params)`` — batched evaluation, returns a
              ``torch.Tensor``.

        solver : str, optional
            ``"batched"`` → :meth:`TorchExpr.torch_integrate_batched`
            (chunked tensor-product quadrature), or ``"vectorized"`` →
            :meth:`TorchExpr.torchquad_integrate_vectorized`.
        change_of_variables_method, cov_eps
            Passed to :meth:`torchify` for each integral term.
        **solver_kwargs
            Forwarded to the chosen solver, e.g. ``N=1001``,
            ``method="gauss-legendre"``, ``chunk_size_points=4096``,
            ``dtype=torch.float64``.

        Returns
        -------
        sympy.Float | torch.Tensor
            Scalar result for scalar parameters, a tensor with the batch shape
            of the parameters otherwise.

        USAGE
        =====
        1. Scalar evaluation::

            ````python
            lt = TorchSymPy()
            result = lt.eval_numeric(expr, solver="batched",
                                     N=1001, method="gauss-legendre")
            ````

        2. Batched evaluation over a mesh::

            ````python
            X, Y = torch.meshgrid(X_t, Y_t, indexing="ij")
            Z = lt.eval_numeric(intensity_eq, params_values=[X, Y],
                                solver="batched", N=1001,
                                method="gauss-legendre", chunk_size_points=4096)
            ````
        """
        if solver not in ("batched", "vectorized"):
            raise ValueError(f"Unknown solver '{solver}'. Expected 'batched' or 'vectorized'.")

        if isinstance(expr, Eq):
            expr = expr.rhs

        free_symbols = sorted(expr.free_symbols, key=lambda symbol: symbol.name)
        params_list = self._split_params_for_lambdify(params_values, len(free_symbols))

        is_batched = any(
            (torch.is_tensor(p) and p.ndim > 0) or (isinstance(p, np.ndarray) and p.ndim > 0)
            for p in params_list
        )

        def solve_integral(integral_node):
            """Evaluate a single Integral.

            A separable integrand is factorised into 1D integrals that are
            refined until they converge (this is what keeps oscillatory
            diffraction integrands correct); a non-separable but shift
            invariant one (Rayleigh-Sommerfeld) is evaluated through a single
            refined cumulative integral of its kernel; everything else goes
            through the original torchify + tensor-product solver path.
            """
            separable_value = self._evaluate_separable_integral(
                integral_node, free_symbols, params_list, **solver_kwargs
            )
            if separable_value is not None:
                return separable_value

            shift_invariant_value = self._evaluate_shift_invariant_integral(
                integral_node, free_symbols, params_list, **solver_kwargs
            )
            if shift_invariant_value is not None:
                return shift_invariant_value

            texpr = self.torchify(
                integral_node,
                params=free_symbols,
                change_of_variables_method=change_of_variables_method,
                cov_eps=cov_eps,
            )
            if solver == "batched":
                re_part, im_part = texpr.torch_integrate_batched(
                    params_values=params_values, **solver_kwargs
                )
                if torch.any(im_part != 0):
                    return re_part + 1j * im_part
                return re_part
            return texpr.torchquad_integrate_vectorized(
                params_values=params_values, **solver_kwargs
            )

        def is_zero_integrand(integral_node) -> bool:
            function = integral_node.function
            return function.is_zero is True or function == 0

        def reraise_if_out_of_memory(exc: BaseException, integral_node) -> None:
            """Never hide a failed allocation behind the zero fallback.

            Returning zeros for an integral that ran out of memory turns a hard
            failure into silently wrong data (a flat image), which is much worse
            than an exception; large parameter meshes therefore get an
            actionable error instead.
            """
            if not _is_out_of_memory_error(exc):
                return
            n_values = max(
                (
                    int(np.prod(np.shape(p)))
                    for p in params_list
                    if torch.is_tensor(p) or isinstance(p, np.ndarray)
                ),
                default=1,
            )
            raise MemoryError(
                f"Ran out of memory while integrating {integral_node} for "
                f"~{n_values} parameter values with solver='{solver}' "
                f"({solver_kwargs}). The integrand is evaluated at "
                "(quadrature points x parameter values) at once, so memory "
                "grows with the product of the two. Fixes, cheapest first: "
                "(a) pass chunk_size_params=<e.g. 4096>, or "
                "chunk_size_params='auto' with solver='vectorized', to cap the "
                "number of parameter values evaluated simultaneously (both "
                "solvers support it; solver='batched' also takes "
                "chunk_size_points), "
                "(b) lower N, (c) use a coarser parameter mesh, or (d) split "
                "the mesh yourself and concatenate the results."
            ) from exc

        if not is_batched:
            # Scalar path: walk the tree, replace integrals by Floats, rebuild.
            def walk_scalar(node):
                if isinstance(node, Integral):
                    if is_zero_integrand(node):
                        return SympyFloat(0.0)
                    try:
                        value = solve_integral(node)
                    except (IndexError, RuntimeError, MemoryError) as exc:
                        reraise_if_out_of_memory(exc, node)
                        logger.warning("integral {} could not be evaluated ({}); using 0", node, exc)
                        value = 0.0
                    if torch.is_tensor(value):
                        if value.numel() != 1:
                            raise RuntimeError(
                                f"Expected a scalar result for {node}; got shape {tuple(value.shape)}"
                            )
                        value = value.item()
                    if isinstance(value, complex):
                        if abs(value.imag) < 1e-12:
                            return SympyFloat(value.real)
                        return SympyFloat(value.real) + sympy_I * SympyFloat(value.imag)
                    return SympyFloat(value)
                if not node.args:
                    return node
                return node.func(*[walk_scalar(arg) for arg in node.args])

            return walk_scalar(expr)

        # Batched path: evaluate the integrals, then lambdify the outer expression.
        dummy_map: Dict[Symbol, Any] = {}

        def walk_batched(node):
            if isinstance(node, Integral):
                if is_zero_integrand(node):
                    return SympyFloat(0.0)
                try:
                    value = solve_integral(node)
                except (IndexError, RuntimeError, MemoryError) as exc:
                    reraise_if_out_of_memory(exc, node)
                    logger.warning("integral {} could not be evaluated ({}); using 0", node, exc)
                    value = torch.zeros(())
                dummy = Symbol(f"_en_{len(dummy_map)}")
                dummy_map[dummy] = value
                return dummy
            if not node.args:
                return node
            return node.func(*[walk_batched(arg) for arg in node.args])

        expr_with_dummies = walk_batched(expr)

        args_symbols = free_symbols + list(dummy_map)
        args_values = [torch.as_tensor(p) for p in params_list] + list(dummy_map.values())

        outer_func = self._lambdify(
            tuple(args_symbols), expr_with_dummies, self._default_modules(), cacheable=True
        )
        result = outer_func(*args_values)

        if torch.is_tensor(result) and torch.is_complex(result):
            if torch.allclose(result.imag, torch.zeros_like(result.imag), atol=1e-12):
                result = result.real

        return result
