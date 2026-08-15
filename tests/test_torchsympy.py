#!/usr/bin/env python
# coding: utf-8

# # libtorch integration test
# 
# This notebook validates that `libtorch.py` integrates SymPy expressions using Torch + torchquad.
# 
# It is intentionally minimal and should run on CPU-only setups; CUDA is optional.

# In[ ]:


import sys
from pathlib import Path
import torch

def _find_repo_root_with_src(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "src" / "torchsympy.py").exists():
            return p
    raise FileNotFoundError("Could not find 'src/libtorch.py' above current working directory")

ROOT = _find_repo_root_with_src(Path.cwd().resolve())
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

print("cwd:", Path.cwd())
print("src:", SRC)
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("cuda:", torch.version.cuda)


# In[2]:


import math
import sympy as sp
import torchsympy
torch.set_default_dtype(torch.float64)
lt = torchsympy.TorchSymPy()
print("sympy:", sp.__version__)


# ### Integral expressions for both tests
# 
# Test 1 (complex Gaussian, over (-∞, ∞))  
# Let x, k ∈ ℝ.  
# $$
# \int_{-\infty}^{\infty} e^{-x^{2}} e^{i k x}\,dx = \sqrt{\pi}\, e^{-k^{2}/4}.
# $$
# 
# Test 2 (semi-infinite exponential, over (0, ∞))  
# Let x ∈ ℝ.  
# $$
# \int_{0}^{\infty} e^{-x}\,dx = 1.
# $$

# In[3]:


# Test 1: complex Gaussian integral over (-oo, oo)
# ∫ exp(-x^2) * exp(I*k*x) dx = sqrt(pi) * exp(-k^2/4)

x, k = sp.symbols("x k", real=True)
integrand = sp.exp(-x**2) * sp.exp(sp.I * k * x)

texpr = lt.torchify(
    integrand,
    variables=[x],
    limits=[(x, -sp.oo, sp.oo)],
    params=[k],
)


# In[4]:


texpr.domain


# In[5]:


texpr.dim


# In[6]:


re, im = texpr.torchquad_integrate(params_values=[1.0], N=151)


# In[7]:


expected = math.sqrt(math.pi) * math.exp(-(1.0**2) / 4.0)
re_f = float(re.detach().cpu())
im_f = float(im.detach().cpu())

print("re:", re_f, "expected:", expected, "abs err:", abs(re_f - expected))
print("im:", im_f)

assert abs(re_f - expected) < 5e-3
assert abs(im_f) < 5e-3


# In[8]:


# Test 2: semi-infinite integral over (0, oo)
# ∫_0^∞ exp(-x) dx = 1

x = sp.Symbol("x", real=True)
expr = sp.exp(-x)

texpr = lt.torchify(
    expr,
    variables=[x],
    limits=[(x, 0, sp.oo)],
    params=[],
)

re, im = texpr.torchquad_integrate(params_values=[], N=151)
re_f = float(re.detach().cpu())
im_f = float(im.detach().cpu())

print("re:", re_f, "expected:", 1.0, "abs err:", abs(re_f - 1.0))
print("im:", im_f)

assert abs(re_f - 1.0) < 5e-3
assert abs(im_f) < 5e-3

