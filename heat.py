import numpy as np
import matplotlib.pylab as p

def dudt(U, A):
    return np.dot(A, U)

def main():
    # u_t = u_xx
    dx = .1
    dt = .5
    timesteps = 100000

    x = np.arange(-10,10,dx)
    m = len(x)
    kappa = 50

    # u''(x) = (u(x + dx) - 2u(x) + u(x - dx)) / dx^2
    ones = lambda x: np.ones(x)
    A = np.diag(ones(m-1),k=-1) + -2*np.diag(ones(m)) + np.diag(ones(m-1),k=1)
    A *= kappa*(dx**2)

    U = 0*ones(m)
    for i in xrange(0,m):
        if x[i] > -2 and x[i] < 2:
            U[i] = 1

    p.ion()
    lines, = p.plot(x,U)
    for n in xrange(0,timesteps):
        U = U + dt*dudt(U,A)
        
        if n % 100 == 0:
            lines.set_ydata(U)
            p.draw()
    p.show()

if __name__ == "__main__":
    main()
