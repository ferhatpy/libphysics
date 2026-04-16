from libphysics.quantum_mechanics import *
from sympy import *

oqmec = quantum_mechanics()
print("qho.a type:", type(oqmec.qho.a))
print("qho.ad type:", type(oqmec.qho.ad))
print("qho.a:", oqmec.qho.a)
print("qho.ad:", oqmec.qho.ad)

n = symbols('n', positive=True, integer=True)
hbar, w = symbols('hbar w', real=True, positive=True)

Hp = S(1)/4*hbar*w*(oqmec.qho.a + oqmec.qho.ad)**4
print("Hp type:", type(Hp))
print("Hp:", Hp)
