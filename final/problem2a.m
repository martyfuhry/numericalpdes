% Problem 2 (a)
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

dx = (b - a) / (m - 1);    % space discretization
dt = 0.9*dx;
x  = [a:dx:b]';      % grid

k1 = [0:m/2-1, 0]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

U = exp(-x.^2);
A = diag(ones(m-1,1),-1) - diag(ones(m-1,1),1);
A(1,m) = 1;
A(m,1) = -1;


for t = 1:timesteps        % loop through the timesteps
    % RK4 time differencing
    Unonlinear = A/(4*dx) * (U.^2);
    k1 = dt * ifft(-i*k.*fft(Unonlinear) + (k.^2 - k.^4).*fft(U));
    Ustar = U + 0.5*k1;
    k2 = dt * ifft((k.^2 - k.^4).*fft(Ustar));
    Ustar = U + 0.5*k2;
    k3 = dt * ifft((k.^2 - k.^4).*fft(Ustar));
    Ustar = U + k3;
    k4 = dt * ifft((k.^2 - k.^4).*fft(Ustar));

    % compute the result for this time step
    U += 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    U
    plot(x,U);             % plot the result
    axis([a,b,0,1])
    drawnow;
endfor
input("Press any key to continue.")
