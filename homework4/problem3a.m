% Problem 3 (a)
% Marty Fuhry
% 3/28/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implement the pseudo-spectral method with leapfrog time differencing for the %
% linear advection equation                                                    %
%  u_t + a u_x = 0.                                                            %
% Use m = 128, a = 1, and u(x,0) = sech(10(x - pi)), choose an appropriate     %
% time step dx for stability. Initialize u_1(x) as you did in Assignment 2.    %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set these pde variables
m = 128;
a = 1;

% let's make the interval [0,2pi]
dx = 2*pi / (m - 1);
x = [0: dx: 2*pi]';

% and set our timesteps
dt = 0.9*dx/(pi * a);
timesteps = 2*pi / dt;

% set the initial condition
u = sech(10.*(x.- pi));

% calculate the first time step
A = diag(ones(m-1,1), 1) - diag(ones(m-1,1), -1);
mu = a * dt / dx;
uprev = u;
u = u - mu*A*u;

% exact second time step
%u = sech(10.*((x - a*dt) - pi));

% for octave's stupid fourier coefficient indexing
k1 = [0:m/2-1, 0]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

for t = 1:timesteps
    % store our previous time step
    utemp = u;

    % leapfrog time differencing
    u = uprev - 2*a*dt*ifft(i*k.*fft(u));

    % use the previous time step next time
    uprev = utemp;
    
    % plot it
    plot(x,u);
    legend("Pseudo-Spectral Method Approximation");
    if t == 1
        title("Problem 3 (a): t = 0");
        print("problem3a1.png");
    elseif t == ceil(timesteps/2)
        title("Problem 3 (a): t = pi");
        print("problem3a2.png");
    elseif t == floor(timesteps)
        title("Problem 3 (a): t = 2*pi");
        print("problem3a3.png");
    endif
    axis([0, 2*pi, -1, 1])
    drawnow;
endfor

input("Press any key to continue.");

