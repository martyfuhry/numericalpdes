% Problem 4 (b)
% Marty Fuhry
% 2/14/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the advection equation                                              %
%  u_t + a u_x = 0                                                             %
% on the interval 0 < x < 1 with a = 1.                                        %
% Use the upwind method                                                        % 
%  u_j^n+1 = u_j^n - mu*(u_j^n - u_j-1^n)                                      %
% with m=128 grid points and the intial condition                              %
%  u(x, 0) = sech^2(20x - 10)                                                  %
% with periodic boundary conditions                                            %
%  u(0, t) = u(1, t).                                                          %
% Demonstrate that this method is convergent of O(dt + dx).                    %
% Pick dt and dx so that                                                       %
%  a * dt / dx = 0.5.                                                          %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Octave is dumb; tell it this is not a function file.
1;

% exact solution is u(x,t) = sech^2(20(x - t) - 10)
function [u] = uexact(x,t)
    u = sech(20*(x - t) - 10).^2;
endfunction
    
% set up some PDE variables
alpha = 1;
a = 0;
b = 1;

dtList = zeros(1,11);
dt1 = 0.05;
for i = 1:11
    dtList(i) = dt1;
    dt1 /= 2;
endfor

% loop through the time values
for dt = [dtList]

    % discretize the grid
    dx = alpha * dt / 0.5;
    x = [a-dx:dx:b+dx]';
    m = size(x)(1) - 2

    timesteps = 10;   % just go for 10 timesteps
    Uexact    = uexact(x, timesteps * dt);
    mu = alpha*dt/dx; % for the FD method

    % define the FD matrix 
    A = sparse(diag(ones(m+2,1)) - diag(ones(m+1,1),-1));
    A(1,m) = -1;    % with periodic 
    A(m,1)  = 1;     % boundary conditions
    U = sech(20*x - 10).^2; % and this initial condition

    % start the FD method
    for t = 1:timesteps
        U += -mu*A*U;
    endfor

    % calculate the norm of the absolute error
    E = norm(Uexact - U, 2)*sqrt(dx);
    printf("dt = %e, dx = %e\n", dt, dx);
    printf("    Error = %e \n\n", E);
endfor
