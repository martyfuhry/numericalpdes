% Problem 3 (a)
% Marty Fuhry
% 2/13/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the heat equation                                                   %
%  u_t = k u_xx                                                                %
% over 0 < x < 1, k = 1                                                        %
% Use Method 1 to solve this over m=128 grid points.                           %
% Forward in time, centered in space:                                          %
%  u_j^n+1 = u_j^n + k dt / dx * ( u_j-1^n - 2 u_j^n + u_j+1^n)                %
% Assume Dirichlet boundary conditions                                         %
%  u(0,t) = u(1,t) = 0                                                         %
% and initial conditions                                                       %
%  u(x,0) = sech^2(200x - 100)                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
kappa = 1;
m = 128;
a = 0;
b = 1;
timesteps = 100;

dx = (b - a) / (m - 1);    % space discretization
x  = [a:dx:b]';            % grid

% tridiagaonal matrix with [1 -2 1] 
A = -2*diag(ones(m,1)) + 1*diag(ones(m-1,1),-1) + 1*diag(ones(m-1,1),1);

% loop through the different time discretizations
for dt = [0.00004, 0.00003, 0.000015]
    U  = sech(200*x - 100).^2; % reset initial conditions
    mu = kappa * dt / dx^2; 
    for t = 1:timesteps   % loop through the timesteps
        U += mu*A*U;      % perform the FD method
        plot(x,U);        % plot the result
        axis([a,b,0,1])
        drawnow;
    endfor
    input("Press any key to continue.")
endfor
