% Problem 2 (c)
% Marty Fuhry
% 4/11/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 128;
a = -20;
b = 20;
timesteps = 1000;

k1 = [0:m/2]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

dx = (b - a) / (m - 1);    % space discretization
kmax = m/2;
dt = 2.5 / (sqrt(kmax^2 + kmax^4));
dt = 2.5 / kmax^2;

x  = [a:dx:b]';      % grid

U = exp(-x.^2);
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
A(1,m) = -1;
A(m,1) = 1;

A *= 1/(4 * dx);

for t = 1:ceil(50/dt)
    % RK4 time differencing
    q1 = dt * (-A * (U.^2))            + dt * ifft((k.^2).*fft(U));
    q2 = dt * (-A * ((U + 0.5*q1).^2)) + dt * ifft((k.^2).*fft(U + 0.5*q1));
    q3 = dt * (-A * ((U + 0.5*q2).^2)) + dt * ifft((k.^2).*fft(U + 0.5*q2));
    q4 = dt * (-A * ((U + 1.0*q3).^2)) + dt * ifft((k.^2).*fft(U + 1.0*q3));

    % compute the result for the first time splitting method
    Ustar = U + 1/6 * (q1 + 2*q2 + 2*q3 + q4);

    % backward euler
    U = ifft(1./(1+dt*k.^4).*fft(Ustar));

    plot(x,U);             % plot the result
    axis([a,b,-1,1]);
    drawnow;
endfor
