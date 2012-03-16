% Problem 2
% Marty Fuhry
% 3/15/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a finite volume code for Burgers equation                              %
%  u_t + u u_x = 0                                                             %
% and compare the performance of Lax-Friedrichs                                %
%                                                                              %
% Lax-Wendroff                                                                 %
%                                                                              %
% and Godunov                                                                  %
%                                                                              %
% schemes for the fluxes.                                                      %
% Solve the Riemann problem (u_l, u_r) =                                       %
%  (0,1), (1,0), (-1,1), and (1,-1).                                           %
% Discuss your results and how the different methods compare (see paper).      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Octave is dumb; tell it this is not a function file
1;

% flux function
function [f] = f(u)
    f = 0.5*u.^2;
endfunction

% set up some PDE variables
m = 50;
a = -1;
b = 1;

% discretize the grid
dx = (b - a) / (m - 1);
x = [a-dx:dx:b+dx]';

% keep the courant number 1 so that there is no dispersion
dt = 0.5*dx;
timesteps = 1/dt; 
global mu = dt/dx;

u_lefts  = [0, 1, -1, 1]; 
u_rights = [1, 0, 1, -1];
for k = 1:4
    u_l = u_lefts(k);
    u_r = u_rights(k);
    u = zeros(m+2,1);

    for j = 1:m+2
        if x(j) <= 0
            u(j) = u_l;
        else
            u(j) = u_r;
        endif
    endfor

    u_lf = u;
    u_lw = u;
    u_g  = u;

    % start the FD method
    for t = 1:ceil(timesteps)
        % inflow conditions
        u_lf(1) = u_l;
        u_lw(1) = u_l;

        %%%%%%%%%%%%%%%%%%%%%%%
        % Lax-Friedrichs
        %%%%%%%%%%%%%%%%%%%%%%%
        %F = [0; f([u_lf(3:m+2)]) - f([u_lf(1:m)]); 0];
        %A = diag(ones(m+1,1),-1) + diag(ones(m+1,1),1);
        %A(m+2, m+1) = 2;

        %u_lf = 0.5*A*u_lf - 0.5*mu*F;

        %%%%%%%%%%%%%%%%%%%%%%%
        % Lax-Wendroff
        %%%%%%%%%%%%%%%%%%%%%%%
        F_left  = [0; f(u_lw(2:m+2)) - f(u_lw(1:m+1))];
        F_right = [f(u_lw(2:m+2)) - f(u_lw(1:m+1)); 0];
        A_right = diag(ones(m+2,1)) + diag(ones(m+1,1),1);
        A_left  = diag(ones(m+2,1)) + diag(ones(m+1,1),-1);
        A_left(1,1)      = 2;
        A_right(m+2,m+2) = 2;
        ustar_left  = 0.5*A_left *u_lw - 0.5*mu*F_left;
        ustar_right = 0.5*A_right*u_lw - 0.5*mu*F_right;

        u_lw = u_lw - mu*(f(ustar_right) - f(ustar_left));
        

        %%%%%%%%%%%%%%%%%%%%%%%
        % Godunov
        %%%%%%%%%%%%%%%%%%%%%%%

        plot(x,u_lw, 'r'); 
        axis([a,b,-1.5,1.5])
        drawnow
        clf;
    endfor
    input("Press any key to continue.");
endfor
