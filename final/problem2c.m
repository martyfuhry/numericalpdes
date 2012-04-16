% Problem 2 (c)
% Marty Fuhry
% 4/11/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use time splitting to derive a method that permits a larger dt and repeat    %
% your simulation from (b).                                                    %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 128;
a = -20;
b = 20;

k1 = [0:m/2]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

dx = (b - a) / (m - 1);    % space discretization
kmax = m/2;
dt = 0.24;

x  = [a:dx:b]';      % grid

U = exp(-x.^2);
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
A(1,m) = -1;
A(m,1) = 1;

A *= 1/(4 * dx);

K = (2*pi/40)^2 * k.^2 - (2*pi/40)^4 * k.^4;

for t = 1:ceil(50/dt)
    % RK4 time differencing
    q1 = -dt * A * (U.^2);
    q2 = -dt * A * ((U + 0.5*q1).^2);
    q3 = -dt * A * ((U + 0.5*q2).^2);
    q4 = -dt * A * ((U + 1.0*q3).^2);

    % compute the result for the first time splitting method
    Ustar = U + 1/6 * (q1 + 2*q2 + 2*q3 + q4);

    % trapezoidal
    U = real(ifft((1 + dt/2*K)./(1 - dt/2*K).*fft(Ustar)));
        plot(x,U);             % plot the result
        axis([a,b,-5,10]);
        title("Problem 2 (c), t = 50");
        legend("Pseudo-Spectral Method Approximation");
        %print("problem2c.png");
        drawnow;
    t*dt

endfor

input("Press any key to continue.");
