# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Usage:
======    

In a jupyter-notebook:
----------------------    
python2.#
execfile('/media/veracrypt1/python/projects/libpython/src/libsympy.py')
python3.#
exec(open('libnotebook.py').read())

In a python module:
-------------------
# Import path for library functions.
import sys
lstPaths = ["/media/veracrypt1/python/projects/libpython/src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
from libnotebook import *
"""
import matplotlib_inline.backend_inline
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import sys
#from libsympy import *

# Initiate rendering Latex, HTML, Math etc.
from IPython import get_ipython
from IPython.display import display, HTML, Latex, Math
from sympy import *
from sympy.interactive import printing
from sympy.vector import *
printing.init_printing()

# Set matplotlib plot preferences.
get_ipython().run_line_magic('matplotlib', 'inline') # %matplotlib inline
get_ipython().run_line_magic('precision', '%.4g')
plt.rcParams["font.family"] = "serif"
plt.rcParams['figure.figsize'] = (5, 3) # Resize figures. # rc('figure', figsize=(5, 4), dpi=100)
plt.rcParams['figure.dpi'] = 100
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "10"
#plt.rcParams["text.usetex"] = True # Render texts of the plots as TeX.
matplotlib_inline.backend_inline.set_matplotlib_formats('svg') # set_matplotlib_formats('svg')
#!mkdir output # Create a directory under content
#plt.style.use('ggplot')

print("libnotebook is loaded.")
