% Problem 4 (d)
% Marty Fuhry
% 2/16/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the advection equation                                              %
%  u_t + a u_x = 0                                                             %
% on the interval 0 < x < 1 with a = 1.                                        %
% Use leapfrog and second-order centered space differencing                    % 
%  u_j^n+1 = u_j^n-1 - a * dt / dx (u_j+1^n - u_j-1^n)                         %
% with m=128 grid points and the intial condition                              %
%  u(x, 0) = 0                                                                 %
% with the inflow boundary condition                                           %
%  u(0, t) = 1 - cos(t).                                                       %
% First proceed naively and assume the right boundary condition is             %
%  u(1, t) = 0.                                                                %
% Discuss what you find (see assignment paper).                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 128;
a = 0;
b = 1;

% discretize the grid
dx = (b - a) / (m - 1);
x = [a-dx:dx:b+dx]';

% loop through the time values
for dt = [0.00781, 0.00391, 0.00195]
    timesteps = 20/dt; 
    mu = alpha*dt/dx; % courant number

    % define the FD matrix 
    A = diag(ones(m+1,1), 1) - diag(ones(m+1,1), -1);

    % initial conditions with inflow
    Uprev = 0*x;
    Ucurrent = 0*x;
    Ucurrent(1) = 1 - cos(dt); % inflow

    % this is the exact initial condition
    %Ucurrent = 

    % start the FD method
    for t = 2:timesteps
        Utemp    = Ucurrent;
        Ucurrent(1) = 1 - cos(dt*t); % inflow
        % absorbing boundary conditions
        Ucurrent(m+2) = Uprev(m+1);
        Ucurrent = Uprev - mu*A*Ucurrent;
        Uprev    = Utemp;
        if mod(t,10) == 0
            plot(x(2:m-1),Ucurrent(2:m-1)); % and plot it
            axis([a,b,-10,10])
            drawnow;
        endif
    endfor
    input("Press any key to continue.");
endfor
