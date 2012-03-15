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

% loop through the time values
dt = 0.05 * dx;
timesteps = 1/dt; % go until t = 1
mu = dt/(2*dx); % for the FD method
kappa = dt^2/(2*dx^2);

% define the FD matrix 
A = diag(ones(m-1,1),1) - diag(ones(m-1,1),-1);
B = diag(ones(m-1,1),1) - 2*diag(ones(m,1)) + diag(ones(m-1,1),-1);

M = [0 1; 1 0];
[P, L] = eig(M);

N = sech(40*(x - 0.5)); % and this initial condition
U = zeros(m,1);

W = P*[U'; N'];

w1 = W(1,:)';
w2 = W(2,:)';

C = -mu*A + kappa*B;

% start the FD method
for t = 1:timesteps
    w1 = w1 - C*w1;
    w2 = w2 + C*w2;

    W = [w1'; w2'];
    U = P' * W;

    u = U(1,:)';
    n = U(2,:)';

    plot(x,n); % and plot it
    axis([a,b,-1,1])
    drawnow;
endfor
input("Press any key to continue.");
