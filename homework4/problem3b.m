% Problem 3 (b)
% Marty Fuhry
% 3/28/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify your code from part (a) to solve the diffusion equation               %
%  u_t = kappa * u_xx.                                                         %
% Let kappa = 0.01 and use the same grid and initial conditions from (a).      %
% Try a few different time steps and describe what happens. Use Von Neumann    %
% stability analysis to explain your results.                                  %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set these pde variables
m = 128;
kappa = 0.01;

% let's make the interval [0,2pi]
dx = 2*pi / (m - 1);
x = [0: dx: 2*pi]';

% and set our timesteps
dt = dx^2 * pi^2 / sqrt(3);
timesteps = 2*pi / dt;

% set the initial condition
u = sech(10.*(x.- pi));

% calculate the first time step
A = diag(ones(m-1,1), 1) - 2*diag(ones(m,1)) + diag(ones(m-1,1), -1);
mu = kappa * dt / dx^2;
uprev = u;
u = u + mu*A*u;

% for octave's stupid fourier coefficient indexing
k1 = [0:m/2-1, 0]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

for t = 1:timesteps
    % store our previous time step
    utemp = u;

    % leapfrog time differencing
    u = uprev + 2*dt*kappa*ifft(-k.^2.*fft(u));

    % use the previous time step next time
    uprev = utemp;
    
    % plot it
    plot(x,u);
    axis([0, 2*pi, -1, 1])
    drawnow;
endfor

input("Press any key to continue.");

