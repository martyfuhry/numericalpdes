%%%%%%%%%%%alternate%%%%%%%%%%%%%%%%%%
% Problem 1 (a)
% Marty Fuhry
% 2/14/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the shallow water equations                                         %
%  u_t + n_x = 0                                                               %
%  n_t + u_x = 0                                                               %
% on the interval 0 < x < 1.                                                   %
% Use Lax Wendroff                                                             % 
%  u_j^n+1 = u_j^n - 1/(2dx) A (u_j+1^n - u_j-1^n)                             %
%            + dt^2/(dx^2) A (u_j+1^n - 2u_j^n + u_j-1^n)                      %
% with m=128 grid points and the intial condition                              %
%  n(x, 0) = sech(40(x - 0.5))                                                 %
%  u(x, 0) = 0                                                                 %
% with boundary condtions                                                      %
%  u(0, t) = 0                                                                 %
%  n(1, t) = 0                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
m = 128;
a = 0;
b = 1;

% discretize the grid
dx = (b - a) / (m - 1);
x = [a:dx:b]';

% initial conditions
u = zeros(m,1);
n = sech(40*(x - 0.5)); 

% loop through the time values
dt = dx;
timesteps = 1/dt; % go until t = 1
mu = dt/dx; 

% define the FD matrices
A = diag(ones(m-1,1),1) - diag(ones(m-1,1),-1);
B = diag(ones(m-1,1),1) - 2*diag(ones(m,1)) + diag(ones(m-1,1),-1);

% this is our FD matrix
%C = -mu/2.*A + mu^2/2.*B;

% start the FD method
for t = 1:timesteps
    %u(1) = 0;
    %n(1) = n(2);
    %u(m) = u(m-1);
    %n(m) = 0;
    % perform the FD method
    u = u + (-mu/2*A + mu^2/2*B)*n;
    n = n + (-mu/2*A + mu^2/2*B)*u;

    %u = unew;
    %n = nnew;

    % plot u and n
    hold on
    plot(x,u, 'r'); 
    plot(x,n, 'b'); 
    axis([a,b,-1,1])
    drawnow;
    clf;
endfor
input("Press any key to continue.");
