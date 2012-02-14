% Problem 3 (d)
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
% Assume boundary conditions                                                   %
%  u_x(0, t) = -1                                                              %
% and                                                                          %
%  u(0,t) = u(1,t) = 0                                                         %
% and initial conditions                                                       %
%  u(x,0) = sech^2(200x - 100)                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Octave is dumb; make sure this isn't a function file.
1;

function [k] = kappa(x)
    k = 0.1 + 0.5*(1 + tanh(6*x - 3));
endfunction

% set up some PDE variables
m = 128;
a = 0;
b = 1;
timesteps = 10000;

dx = (b - a) / (m - 1); % space discretization
x  = [a:dx:b]';         % grid

% tridiagaonal matrix with [1 -2 1] 
A = -2*diag(ones(m,1)) + 1*diag(ones(m-1,1),-1) + 1*diag(ones(m-1,1),1);
A = sparse(A);

% approximate kappa derivative and put into matrix
alpha = kappa(x(1:m-1).+dx) - kappa(x(1:m-1).-dx);
B = 1/4*(diag(alpha, 1) - diag(alpha, -1));

% multiply the j'th row of A with kappa(x_j)
for j = 1:m
    A(j,:)*=kappa(x(j));
endfor

I = diag(ones(m,1));

% loop through the different time discretizations
for dt = [0.0004]
    C = dt/dx^2 * (A + B);
    U  = sech(200*x - 100).^2; % reset initial conditions
    for t = 1:timesteps        % loop through the timesteps
        U = (I + C)*U;         % perform the FD method
        U = (I - C)\U;
        plot(x,U);             % plot the result
        axis([a,b,0,.1])
        drawnow;
    endfor
    input("Press any key to continue.")
endfor
