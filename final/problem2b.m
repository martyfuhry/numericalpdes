% Problem 2 (b)
% Marty Fuhry
% 4/11/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the problem                                                            %
%  u_t + u u_x = -u_xx - u_xxxx                                                %
% over the domain [-20,20] with Runge-Kutta 4 in time and a pseudo-spectral    %
% space discretization.                                                        %
% Discretize the nonlinear term by approximating the x-derivative of u^2/2.    %
% Use m = 128 and choose dt as large as possible.                              %
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
kmax = 2*pi/40 * m/2;
dt = 2.5 / (sqrt(kmax^2 + kmax^4 - 2*kmax^6 + kmax^8));
x  = [a:dx:b]';      % grid

U = exp(-x.^2);
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
A(1,m) = -1;
A(m,1) = 1;

A *= 1/(4 * dx);
A = sparse(A);

dt
K = (2*pi/40)^2 * k.^2 - (2*pi/40)^4 * k.^4;

for t = 1:ceil(50/dt)

    % RK4 time differencing
    q1 = dt * (-A*(U.^2) + real(ifft(K.*fft(U))));
    q2 = dt * (-A*((U + 0.5*q1).^2) + real(ifft(K.*fft(U + 0.5*q1))));
    q3 = dt * (-A*((U + 0.5*q2).^2) + real(ifft(K.*fft(U + 0.5*q2))));
    q4 = dt * (-A*((U + 1.0*q3).^2) + real(ifft(K.*fft(U + 1.0*q3))));

    % compute the result for this time step
    U += 1/6 * (q1 + 2*q2 + 2*q3 + q4);
    if (mod(t,100) == 0)
        plot(x,U);             % plot the result
        t*dt
        title("Problem 2 (b), t = 50");
        legend("Pseudo-Spectral Method Approximation");
        axis([a,b,-5,10]);
        drawnow;
    endif
endfor

plot(x,U);             % plot the result
title("Problem 2 (b), t = 50");
legend("Pseudo-Spectral Method Approximation");
axis([a,b,-5,10]);
print("problem2b.png");
drawnow;


