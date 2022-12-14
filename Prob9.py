import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt


# Import double-slit, single-slit and triple-slit
U2 = pa.cx_cube()
U2.load("UsForProb8.bin")
U1 = pa.cx_cube()
U1.load("UsProb9Single.bin")
U3 = pa.cx_cube()
U3.load("UsProb9Triple.bin")

K, J, I = np.shape(U1)
p1 = np.zeros( (I,J,K), dtype=complex )
p2 = np.zeros( (I,J,K), dtype=complex )
p3 = np.zeros( (I,J,K), dtype=complex )
for k in range(K):
    for j in range(J):
        for i in range(I):
            p1[i,j,k] = np.conj(U1[i,j,k])*U1[i,j,k] # Probability
            p2[i,j,k] = np.conj(U2[i,j,k])*U2[i,j,k]
            p3[i,j,k] = np.conj(U3[i,j,k])*U3[i,j,k]
p1 = p1.real # imaginary part is zero, but this keeps pyplot from complaining.
p2 = p2.real
p3 = p3.real
Dt = 2.5e-5 # Time-step
T = 0.002 # Time-period
t = np.arange(0, T, Dt)


h = 0.005 # Length-step
x = np.arange(h,1-h, h)
y = np.arange(h,1-h, h)

# Normalize 
p1_normed = p1[160-1,:,-1]/sum(p1[160-1,:,-1])
p2_normed = p2[160-1,:,-1]/sum(p2[160-1,:,-1])
p3_normed = p3[160-1,:,-1]/sum(p3[160-1,:,-1])

plt.figure()
plt.subplot(221)
plt.plot(y, p1_normed)
plt.xlabel('y')
plt.ylabel('p(y)')
plt.subplot(222)
plt.plot(y, p2_normed)
plt.xlabel('y')
plt.ylabel('p(y)')
plt.subplot(223)
plt.plot(y, p3_normed)
plt.xlabel('y')
plt.ylabel('p(y)')
plt.show()