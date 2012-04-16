%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GALERKIN SPECTRAL METHOD TEST                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 128;

a = -20;
b = 20;

% the grid
dx = (b - a) / (m - 1); 
x  = [a:dx:b]';      

% fourier modes
k = [-m/2+1:m/2]';

% time step
kmax = m/2;
dt = 0.01*dx;

% initial condition
U = exp(-x.^2);

% need to get coefficients from ICs
for j = 1:m
    for l = 1:m
        A(j,l) = exp(i*2*pi/40*k(l)*x(j));
    endfor
endfor

% solve the system to get the coefficients
C = A\U;
    
U = 0;
for j = 1:m
    U += C(j)*exp(i*2*pi/40*k(j)*x);
endfor

% for the nonlinear term
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
% and boundary conditions
A(1,m) = -1;
A(m,1) = 1;

A *= 1/(4 * dx);
A = sparse(A);

for t = 1:ceil(50/dt)

    C = C - 2*pi/40*dt*(i*k - k.^2).*C;

    U = 0;
    for j = 1:m
        U += C(j)*exp(i*2*pi/40*k(j)*x);
    endfor

    t*dt
    plot(x,U);             % plot the result
    axis([a,b,-1,1]);
    drawnow;
endfor
