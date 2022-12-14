import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

U0 = pa.cx_cube() # Import wavefunction, no barrier. 
U0.load("Usv0.bin")

K, J, I = np.shape(U0)
P0 = np.zeros(K, dtype=complex)
for k in range(K):
    for j in range(J):
        for i in range(I):
            P0[k] += np.conj(U0[i,j,k])*U0[i,j,k] # Sum of probability-function 
                                                    # at time k. Should be 1. 
Dt = 2.5e-5 # Time step
T = 0.008 # Total period 
t = np.arange(0, T, Dt) 

plt.figure()
plt.plot(t,1-P0.real)
plt.xlabel('t')
plt.ylabel('1 - P(t)')
plt.show()

U1 = pa.cx_cube() # With double-slit. 
U1.load("Usv1e10.bin")

K, J, I = np.shape(U1)
P1 = np.zeros(K, dtype=complex)
for k in range(K):
    for j in range(J):
        for i in range(I):
            P1[k] += np.conj(U1[i,j,k])*U1[i,j,k]


plt.figure()
plt.plot(t,1-P1.real)
plt.xlabel('t')
plt.ylabel('1 - P(t)')
plt.show()