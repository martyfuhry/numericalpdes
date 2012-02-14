% Problem 3 (c)
% Marty Fuhry
% 2/13/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the heat equation                                                   %
%  u_t = k u_xx                                                                %
% over 0 < x < 1, k = 1                                                        %
% Use Crank-Nicolson to solve this over m=128 grid points.                     %
%  (I - mu) U^n+1 = (I + mu) U^n + b
% Assume boundary conditions                                                   %
%  u_x(0, t) = -1                                                              %
% and                                                                          %
%  u(0,t) = u(1,t) = 0                                                         %
% and initial conditions                                                       %
%  u(x,0) = sech^2(200x - 100)                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
kappa = 1;
m = 128;
a = 0;
b = 1;
timesteps = 100000; % this is superfluous

dx = (b - a) / m;    % space discretization
x  = [a:dx:b]';      % grid

% tridiagaonal matrix with [1 -2 1] 
A = -2*diag(ones(m+1,1)) + 1*diag(ones(m,1),-1) + 1*diag(ones(m,1),1);
A = sparse(A);
A(1,:) = [-2, 2, zeros(1,m-1)]; % need the first row to conform to Neumann boundary conditions

B = zeros(m+1, 1); % boundary conditions
I = diag(ones(m+1,1));

% loop through the different time discretizations
% since the first choice is unstable, we'll only look at the third
for dt = [0.00003]
    U  = sech(200*x - 100).^2; % reset initial conditions
    mu = kappa * dt / (2 * dx^2);
    B(1) = mu * 4*dx;          % boundary conditions
    for t = 1:timesteps        % loop through the timesteps
        U = (I + A.*mu)*U + B;       % perform the FD method
        U = (I - A.*mu)\U;
        plot(x,U);             % plot the result
        axis([a,b,0,1])
        % print out and draw every ten timesteps
        if mod(t,10) == 0
            printf("t = %f\n", t*dt);
            drawnow;
        endif
    endfor
    input("Press any key to continue.")
endfor
