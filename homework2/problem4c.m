% Problem 4 (c)
% Marty Fuhry
% 2/14/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the advection equation                                              %
%  u_t + a u_x = 0                                                             %
% on the interval 0 < x < 1 with a = 1.                                        %
% Use leapfrog and second-order centered space differencing                    % 
%  u_j^n+1 = u_j^n-1 - a * dt / dx (u_j+1^n - u_j-1^n)                         %
% with m=128 grid points and the intial conditions                             %
%  u(x, 0) = sech^2(20x - 10)                                                  %
%  and the narrower                                                            %
%  u(x, 0) = sech^2(40x - 20)                                                  %
% with periodic boundary conditions                                            %
%  u(0, t) = u(1, t).                                                          %
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
for dt = [0.0156, 0.00781, 0.00391, 0.00195]
    timesteps = 1/dt; % make it loop once
    mu = alpha*dt/dx; % for the FD method

    % define the FD matrix 
    A = diag(ones(m+1,1), 1) - diag(ones(m+1,1), -1);
    A(1,m+1) = -1;    % with periodic 
    A(m+2,2) = 1;     % boundary conditions

    Uprev = sech(40*x - 20).^2; % and this initial condition first

    % use explicit trapezoidal method for one step
    Ustar = Uprev - mu * A * Uprev;
    Ucurrent = Uprev - 0.5*mu*A*(Ustar + Uprev);

    % start the FD method
    for t = 0:timesteps
        Utemp    = Ucurrent;
        Ucurrent = Uprev - mu*A*Ucurrent;
        Uprev    = Utemp;
        plot(x(2:m-1),Ucurrent(2:m-1)); % and plot it
        axis([a,b,-1,1])
        drawnow;
    endfor
    input("Press any key to continue.");

endfor

