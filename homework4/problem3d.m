% Problem 3 (d)
% Marty Fuhry
% 3/28/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combine your code from parts (a) and (c) to build an FFT method for the      %
% advection-diffusion equation                                                 %
%  u_t + a * u_x = kappa * u_xx.                                               %
% Solve this using a = 1 and kappa = 0.01 and the same grid and initial        %
% condition as in part (a).                                                    %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set these pde variables
m = 128;
kappa = 0.01;
alpha = 1;

% let's make the interval [0,2pi]
dx = 2*pi / (m - 1);
x = [0: dx: 2*pi]';

% and set our timesteps
dt = dx;
timesteps = 2*pi / dt;

% set the initial condition
u = sech(10.*(x.- pi));

% store previous time step 
uprev   = u;

% calculate the first time step
A = diag(ones(m-1,1), 1) - 2*diag(ones(m,1)) + diag(ones(m-1,1), -1);
B = diag(ones(m-1,1), 1) - diag(ones(m-1,1),-1);

u = u + (-alpha*dt/dx*A + kappa*dt/dx^2*B)*u;

% for octave's stupid fourier coefficient indexing
k1 = [0:m/2]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

for t = 1:timesteps
    % store our previous time step
    utemp = u;

    % crank-nicolson time differencing
    u = ifft(((1. - i*k*dt*alpha - kappa*dt*k.^2.))./(1 + i*k*dt*alpha + kappa*k.^2*dt).*fft(uprev));

    % use the previous time step next time
    uprev = utemp;
    
    % plot it
    plot(x,u);
    axis([0, 2*pi, -1, 1])
    drawnow;
    legend("Pseudo-Spectral Method Approximation");
    if t == 1
        title("Problem 3 (d): t = 0");
        print("problem3d1.png");
    elseif t == ceil(timesteps/2)
        title("Problem 3 (d): t = pi");
        print("problem3d2.png");
    elseif t == floor(timesteps)
        title("Problem 3 (d): t = 2*pi");
        print("problem3d3.png");
    endif
 
endfor

input("Press any key to continue.");

