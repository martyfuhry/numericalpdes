% Problem 2 (d)
% Marty Fuhry
% 4/11/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify your code from (b) to implement the Galerkin spectral method.         %
%                                                                              %
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
kmax = 2*pi/40 * m/2;
dt = 2.5 / (sqrt(kmax^2 + kmax^4 - 2*kmax^6 + kmax^8));

% fourier modes
k1 = [0:m/2]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];


% initial condition
U = exp(-x.^2);

% need to get coefficients from ICs
for j = 1:m
    for l = 1:m
        B(j,l) = exp(i*2*pi/40*k(l)*x(j));
    endfor
endfor

% solve the system to get the coefficients
C = B\U;

% for the nonlinear term
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
% and periodic boundary conditions
A(1,m) = -1;
A(m,1) = 1;

A *= 1/(4 * dx);
A = sparse(A);

K = (2*pi/40)^2 * k.^2 - (2*pi/40)^4 * k.^4;

M = 3/2*m + 2;
K1 = [0:M/2]';
K2 = [-M/2 + 1: -1]';
otherK = [K1; K2];
index = find(otherK>-m/2 & otherK<m/2+1);

for t = 1:ceil(50/dt)

    % need to do the convolution for the nonlinear term
    Ck = zeros(M,1);
    Ck(index) = C;
    UUx = ifft(fft((2*pi/40)*i*Ck.*otherK).*fft(Ck));
    UUx = UUx(index);

    % linear diffusion terms
    q1 = -dt*UUx + dt * K.*C;
    q2 = -dt*UUx + dt * K.*(C + 0.5*q1);
    q3 = -dt*UUx + dt * K.*(C + 0.5*q2);
    q4 = -dt*UUx + dt * K.*(C + 1.0*q3);

    % compute the result for this time step
    C += 1/6 * (q1 + 2*q2 + 2*q3 + q4);

    % evaulate U using the coefficients
    U = B*C;

    if (mod(t,1000) == 0)
        t*dt
        plot(x,U);             % plot the result
        axis([a,b,-5,10]);
        title("Problem 2 (d), t = 50");
        legend("Galerkin-Spectral Method Approximation");
        drawnow
    endif
endfor

plot(x,U);             % plot the result
axis([a,b,-5,10]);
title("Problem 2 (d), t = 50");
legend("Galerkin-Spectral Method Approximation");
print("problem2d.png");
drawnow;

