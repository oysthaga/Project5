import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt

# Import potential
V = pa.mat()
V.load("Mat1.bin")

# Define grid
x = np.linspace(0,1,np.shape(V)[0]+1)
y = np.linspace(0,1,np.shape(V)[1]+1)

# Plot
plt.figure()
plt.pcolormesh(x,y, V.t()) # V.t() transposes matrix, because plt.pcolormesh
                            # puts x along vertical axis. 
plt.xlabel('x'); plt.ylabel('y')
plt.show()