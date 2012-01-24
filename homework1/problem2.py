#u_t = -a u_x - B u_xxx

import numpy as np
import matplotlib.pylab as p

h = 0.1
dx = h
dt = .005
x = np.arange(-1-h,1+3*h,h)
m = len(x)
n = 20

a = 2
beta = .1

U = np.empty((n, m))

c = .5
U[0,:] = -2*np.sin(c*x)

A = -1*np.diag(np.ones(m)) + np.diag(np.ones(m-1), k=1)
A[m-1,1] = 1

B = np.diag(-1*np.ones(m-1), k=-1) + 3*np.diag(np.ones(m)) + -3*np.diag(np.ones(m-1), k=1) + 1*np.diag(np.ones(m-2),k=2)
B[0,m-2] = -1
B[m-2, 1] = 1
B[m-1, 2] = 1
B[m-1, 1] = -3


A*=a/dx
B*=beta/(dx**3)

for j in xrange(0,n-1):
    fig = p.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim((-2,2))
    ax.set_autoscaley_on(False)
    p.title("$t="+str(j)+"$")
    q1 = h*(dt*np.dot((A + B), U[j,:]))
    p1 = U[j,:] + q1
    q2 = h*(dt*np.dot((A + B), p1))
    U[j+1,:] = p1 + q2 / 2
    #U[j+1,:] = U[j,:] - dt*np.dot((-A - B), U[j,:])
    ax.plot(x,U[j+1,:])
    p.savefig("fig"+str(j)+".png")
