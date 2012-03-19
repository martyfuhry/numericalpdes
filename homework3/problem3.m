% Problem 2
% Marty Fuhry
% 3/15/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a finite volume code for Burgers equation                              %
%  u_t + u u_x = 0                                                             %
% and compare the performance of Lax-Friedrichs                                %
% Lax-Wendroff                                                                 %
% and Godunov                                                                  %
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

function [m] = minf(u_l, u_r)
    tol = (u_r - u_l) / 10.;
    if abs(tol) <= 10e-16;
        m = f(u_l);
    else
        u = [u_l:tol:u_r];
        m = min(f(u));
    endif
endfunction

function [m] = maxf(u_l, u_r)
    tol = (u_r - u_l) / 10.;
    if abs(tol) <= 10e-16;
        m = f(u_l);
    else
        u = [u_l:tol:u_r];
        m = max(f(u));
    endif
endfunction

% set up some PDE variables
m = 100;
a = -1;
b = 1;

% discretize the grid
dx = (b - a) / (m - 1);
x = [a:dx:b]';
dt = dx;
timesteps = 1/dt; 
mu = dt/dx;

u_lefts  = [0, 1, -1, 1]; 
u_rights = [1, 0, 1, -1];
for k = 1:4
    u_l = u_lefts(k);
    u_r = u_rights(k);
    u = zeros(m,1);

    for j = 1:m
        if x(j) <= 0
            u(j) = u_l;
        else
            u(j) = u_r;
        endif
    endfor

    u_lf = u;
    u_lw = u;
    u_g  = u;

    % start the FV method
    for t = 1:ceil(timesteps)
        % inflow conditions
        u_lf(1) = u_l;
        u_lw(1) = u_l;
        u_g(1)  = u_l;

        %%%%%%%%%%%%%%%%%%%%%%%
        % Lax-Friedrichs      %
        %%%%%%%%%%%%%%%%%%%%%%%
        % calculate the flux
        flux = f(u_lf);
        F = diag(ones(m-1,1),1) - diag(ones(m-1,1),-1);
        F(1,1) = -1;
        F(m,m) = 1;
        A = diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
        A(1, 1)   = 1; % this should be fine
        A(m, m-1) = 2; % need to get outflow right

        u_lf = 0.5*A*u_lf - 0.5*mu*F*flux;

        %%%%%%%%%%%%%%%%%%%%%%%
        % Lax-Wendroff        %
        %%%%%%%%%%%%%%%%%%%%%%%
        % calculate the flux
        flux = f(u_lw);
        F_left  = diag(ones(m,1)) - diag(ones(m-1,1),-1);
        F_right = diag(ones(m-1,1),1) - diag(ones(m,1));
        F_left(1,1)  = 0; % already have inflow 
        F_right(m,m) = 0; % this should be fine
        A_right = diag(ones(m,1)) + diag(ones(m-1,1),1);
        A_left  = diag(ones(m,1)) + diag(ones(m-1,1),-1);
        A_left(1,1)  = 2; % already have inflow
        A_right(m,m) = 2; % need to get outflow right
        ustar_left  = 0.5*A_left *u_lw - 0.5*mu*F_left*flux;
        ustar_right = 0.5*A_right*u_lw - 0.5*mu*F_right*flux;

        u_lw = u_lw - mu*(f(ustar_right) - f(ustar_left));
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Godunov             %
        %%%%%%%%%%%%%%%%%%%%%%%
        F = zeros(m,1);
        F(1) = f(u_l); % inflow
        % calculate flux
        for i = 2:m
            if u_g(i-1) <= u_g(i)
                F(i) = minf(u_g(i-1), u_g(i));   
            else
                F(i) = maxf(u_g(i-1), u_g(i));
            endif
         endfor
        u_g = u_g - mu*([F(2:m); F(m)] - F);

        hold on
        plot(x,u_lf, 'b.-');
        plot(x,u_lw, 'r-'); 
        plot(x,u_g, 'g--'); 
        legend("Lax-Friedrich", "Lax-Wendroff", "Godunov")
        axis([a,b,-1.5,1.5])
        drawnow
        if t == ceil(timesteps/2)
            if k == 1
                print "-S640,480" prob3-k1.png
            elseif k == 2
                print "-S640,480" prob3-k2.png
            elseif k == 3
                print "-S640,480" prob3-k3.png
            elseif k == 4
                print "-S640,480" prob3-k4.png
            endif
        endif
        clf;
    endfor
    input("Press any key to continue.");
endfor
