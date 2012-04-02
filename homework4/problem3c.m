% Problem 3 (c)
% Marty Fuhry
% 3/28/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify your code from part (b) to use Crank-Nicolson                         %
%  u_j^n+1 = u_n^n-1 + kappa * dt * ((u_xx)_j^n-1 + (u_xx)_j^n+1)              %
% Implement this method and use it to re-do your experiments from part (b).    %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set these pde variables
m = 128;
kappa = 0.01;

% let's make the interval [0,2pi]
dx = 2*pi / (m - 1);
x = [0: dx: 2*pi]';

% and set our timesteps
dt = 0.1*dx;
timesteps = 2*pi / dt;

% set the initial condition
u = sech(10.*(x.- pi));

% calculate the first time step
A = diag(ones(m-1,1), 1) - 2*diag(ones(m,1)) + diag(ones(m-1,1), -1);
mu = kappa * dt / dx^2;

% store previous time step variables
uprev = u;

u = u + mu*A*u;

% for octave's stupid fourier coefficient indexing
k1 = [0:m/2-1,0]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

for t = 1:timesteps
    % store our previous time step
    utemp = u;

    % crank-nicolson time differencing
    %u = uprev + kappa * dt * (uxxprev + uxxH
    u = ifft((1. - kappa*dt*k.^2.).*fft(uprev)./(1 + kappa*k.^2*dt));

    % use the previous time step next time
    uprev = utemp;
    
    % plot it
    plot(x,u);
    axis([0, 2*pi, -1, 1])
    legend("Pseudo-Spectral Method Approximation");
    drawnow;
    if t == 1
        title("Problem 3 (c): t = 0");
        print("problem3c1.png");
    elseif t == ceil(timesteps/2)
        title("Problem 3 (c): t = pi");
        print("problem3c2.png");
    elseif t == floor(timesteps)
        title("Problem 3 (c): t = 2*pi");
        print("problem3c3.png");
    endif
endfor

input("Press any key to continue.");

