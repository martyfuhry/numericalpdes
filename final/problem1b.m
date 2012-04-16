% Problem 1 (b)
% Marty Fuhry
% 4/11/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 1000;
a = 0;
b = 1;
timesteps = 1000;

dx = (b - a) / (m - 1);    % space discretization
dt = 1.0*dx;
x  = [a:dx:b]';      % grid

A = 3*diag(ones(m,1)) - 4*diag(ones(m-1,1),-1) + diag(ones(m-2,1),-2);
B = 1*diag(ones(m,1)) - 2*diag(ones(m-1,1),-1) + diag(ones(m-2,1),-2);

mu = alpha * dt / dx;

C = eye(m,m) + -mu/2 * A + mu^2/2 * B;

U = sech(10*(x - 0.5));

for t = 1:timesteps        % loop through the timesteps
    U = C*U;
    plot(x,U);             % plot the result
    axis([a,b,0,1])
    drawnow;
endfor
input("Press any key to continue.")
