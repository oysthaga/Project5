import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
'''
V = pa.mat()
V.load("Mat1.bin")
plt.pcolormesh(V.t())
'''

U1 = pa.cx_cube()
U1.load("Usv1e10.bin")

K, J, I = np.shape(U1)
P = np.zeros(K, dtype=complex)
for k in range(K):
    for j in range(J):
        for i in range(I):
            P[k] += np.conj(U1[i,j,k])*U1[i,j,k]
            
Dt = 2.5e-5 
T = 0.008
t = np.arange(0, T, Dt)
plt.figure()
plt.plot(t,1-P.real)
plt.xlabel('t [s]')
plt.ylabel('1 - P(t)')
plt.show()