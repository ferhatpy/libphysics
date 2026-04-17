import sys
import os

# Ensure libphysics is in path
sys.path.append(os.getcwd())

from libphysics.quantum_optics import oqopt, alpha, beta, k

print("Successfully imported oqopt and symbols.")
print(f"alpha: {alpha}")
print(f"k: {k}")
