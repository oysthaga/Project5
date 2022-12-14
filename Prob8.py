import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

U1 = pa.cx_cube()
U1.load("UsForProb8.bin")

K, J, I = np.shape(U1)
p = np.zeros( (I,J,K), dtype=complex )
for k in range(K):
    for j in range(J):
        for i in range(I):
            p[i,j,k] = np.conj(U1[i,j,k])*U1[i,j,k] # Probability at time and position. 
p = p.real # Discard imaginary part (is zero, but this keeps Pyplot from complaining)
Dt = 2.5e-5 # time-step
T = 0.002 # total period
t = np.arange(0, T, Dt)

# Indices of time-innstances. 
n1 = 0
n2 = 40-1
n3 = 80-1

h = 0.005 # length-step. 
x = np.arange(h,1-h, h)
y = np.arange(h,1-h, h)

plt.figure()
plt.subplot(311)
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.pcolormesh(x,y,p[:,:,n1].transpose())
plt.subplot(312)
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.pcolormesh(x,y,p[:,:,n2].transpose())
plt.subplot(313)
plt.pcolormesh(x,y,p[:,:,n3].transpose())
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.show()

plt.figure()
plt.subplot(311)
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.pcolormesh(x,y,np.real(U1)[n1,:,:].transpose())
plt.subplot(312)
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.pcolormesh(x,y,np.real(U1)[n2,:,:].transpose())
plt.subplot(313)
plt.pcolormesh(x,y,np.real(U1)[n3,:,:].transpose())
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.show()

plt.figure()
plt.subplot(311)
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.pcolormesh(x,y,np.imag(U1)[n1,:,:].transpose())
plt.subplot(312)
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.pcolormesh(x,y,np.imag(U1)[n2,:,:].transpose())
plt.subplot(313)
plt.pcolormesh(x,y,np.imag(U1)[n3,:,:].transpose())
plt.axis('equal')
plt.xlabel('x'); plt.ylabel('y')
plt.show()