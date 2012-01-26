#u_t = -a u_x - B u_xxx

import numpy as np
import matplotlib.pylab as p
import time

def dudt(U, C):
    #dudt = np.concatenate((U[len(U)-2:], U, U[0:1]))
    dudt = np.dot(C, U)
    #return dudt[2:len(U)+2]
    return dudt

dx = .1
dt = dx**3
x = np.arange(-20,20+dx,dx)
m = len(x)
n = int(10/dt)

a = 2
beta = 1

U = np.empty((n, m))

U = 20*np.exp(-x**2)

A = -1*np.diag(np.ones(m-1), k=-1) + 1*np.diag(np.ones(m-1), k=1)
A[0,m-2] = -1
A[m-1, 1] = 1

# u''' = u_j+2 - 2u_j+1 + 2uj-1 -u_j-2 / (2dx^3)
B = -np.diag(np.ones(m-2), k=-2) \
  +  2*np.diag(np.ones(m-1), k=-1) \
  + -2*np.diag(np.ones(m-1), k=1)  \
  +  np.diag(np.ones(m-2),k=2)

B[0,m-2] = 2
B[0,m-3] = -1
B[1,m-2] = -1
B[m-2,1] = 1
B[m-1,2] = 1
B[m-1,1] = -2

A *= a/(2*dx)
B *= beta/(2*dx**3)

C = A + B

p.ion()
fig = p.figure()
ax = fig.add_subplot(111)
ax.set_autoscaley_on(False)
ax.set_ylim((-20,20))
lines, = p.plot(x,U)

for j in xrange(0,n-1):
    k1 = dt * dudt(U, C)
    k2 = dt * dudt(U+k1/2, C)
    k3 = dt * dudt(U+k2/2, C)
    k4 = dt * dudt(U+k3, C)

    U = U + k1/6 + k2/3 + k3/3 + k4/6
    lines.set_ydata(U)
    p.draw()

p.show()
time.sleep(2)



