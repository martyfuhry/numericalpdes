import time
import numpy as np
import matplotlib.pylab as p

# u_t + u u_x = 0

def dudt(U,A):
    return np.dot(A,U)

def main():
    dx = 0.01
    dt = 0.0001
    x  = np.arange(-2,2,dx)
    m  = len(x)
    timesteps = int(1./dt)

    A = np.diag(np.ones(m-1), k=-1) + -1*np.diag(np.ones(m-1),k=1)
    A[0,m-2]  = 1
    A[m-1, 1] = -1

    A *= 1/(2*dx)

    U = 5*np.exp(-x**2/2)

    p.ion()
    lines, = p.plot(x,U)
    for n in xrange(0,timesteps):
        U = U + dt*U*dudt(U,A)
        lines.set_ydata(U)
        p.draw()

main()
