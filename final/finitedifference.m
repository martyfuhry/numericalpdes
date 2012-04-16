%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINITE DIFFERENCE TEST OF PROBLEM 2                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 128;

a = -20;
b = 20;

% the grid
dx = (b - a) / (m - 1); 
x  = [a:dx:b]';      

% time step
kmax = m/2;
dt = 0.001*dx;

% fourier modes
k = [-m/2+1:m/2]';

% initial condition
U = exp(-x.^2);

% for the nonlinear term
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
% and boundary conditions
A(1,m) = -1;
A(m,1) = 1;

A *= -1/(4 * dx);

B = -diag(ones(m-1,1),1) + 2*diag(ones(m,1)) - diag(ones(m-1,1),-1);
B(1,m) = -1;
B(m,1) = -1;

B *= 1/(2*dx^2);

C = -diag(ones(m-2,1),2) + 4*diag(ones(m-1,1),1) - 6*diag(ones(m,1)) + 4*diag(ones(m-1,1),-1) - diag(ones(m-2,1),-2);

C *= 1/(dx^4);

D = sparse(B + C);

for t = 1:ceil(50/dt)

    q1 = dt*A*(U.^2) + dt * D*U;
    q2 = dt*A*((U + 0.5*q1).^2) + dt * D*(U + 0.5*q1);
    q3 = dt*A*((U + 0.5*q2).^2) + dt * D*(U + 0.5*q2);
    q4 = dt*A*((U + 1.0*q3).^2) + dt * D*(U + 1.0*q3);

    % compute the result for this time step
    U += 1/6 * (q1 + 2*q2 + 2*q3 + q4);
    dt*t

    if (mod(t,100) == 0)
        plot(x,U);             % plot the result
        axis([a,b,-1,1]);
        drawnow;
    endif
endfor
