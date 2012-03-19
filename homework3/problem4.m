% Problem 4
% Marty Fuhry
% 3/17/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a finite volume code for linear advection with constant a=1            %
%  u_t + a u_x = 0                                                             %
% and compare the performance of upwinding and Lax-Wendroff                    %
%  F = au_j-1^n + a/2 (1-mu)(u_j^n - u_j-1^n)                                  %
% using flux limiters minmod                                                   %
%  max(0, min(1, r))                                                           %
% superbee                                                                     %
%  max(0, min(1, 2r), min(2, r))                                               %
% and Van Leer                                                                 %
%  (r + abs(r))/(1 + abs(r))                                                   %
% How do these compare (see paper)?                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Octave is dumb; tell it this is not a function file
1;

function [C] = minmod(r)
    m = size(r)(1);
    C = zeros(m,1);
    for k = 1:m
        C(k) = max(0, min(1,r(k)));
    endfor
endfunction

function [C] = superbee(r)
    m = size(r)(1);
    C = zeros(m,1);
    for k = 1:m
        C(k) = max([0, min(1,2*r(k)), min(2,r(k))]);
    endfor
endfunction

function [C] = vanleer(r)
    C = (r + abs(r))./(1 + abs(r));
endfunction

% set up some PDE variables
m = 100;
a = -1;
b = 10;

% discretize the grid
dx = (b - a) / (m - 1);
x = [a:dx:b]';
dt = 0.5*dx;
timesteps = (b-a)/dt; 
mu = dt/dx;

% initial conditions
u = zeros(m,1);
u_l = 2;
u_r = 0;
for j = 1:m
    if x(j) <= 0
        u(j) = u_l;
    else
        u(j) = u_r;
    endif
endfor

% set each seperate methods initial conditions
uminmod   = u;
usuperbee = u;
uvanleer  = u;
uupwind   = u;
uexact    = u;

% set left flux matrices
A_left = diag(ones(m,1)) - diag(ones(m-1,1),-1);
A_left(1,1) = 0;
B_left = diag(ones(m-1,1), -1);
B_left(1,1) = 1;

% set right flux matrices
A_right = diag(ones(m-1,1),1) - diag(ones(m,1));
A_right(m,m) = 0;
B_right = diag(ones(m,1));

A = diag(ones(m-1,1),-1) - diag(ones(m-2,1),-2);
B = diag(ones(m,1)) - diag(ones(m-1,1),-1);

% start the FV method
for t = 1:timesteps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % minmod                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        % get ratios          %
        %%%%%%%%%%%%%%%%%%%%%%%
        r1 = A*uminmod;
        r2 = B*uminmod;
        r2(r2 == 0) = 1e-10;
        r = r1./r2;

        %%%%%%%%%%%%%%%%%%%%%%%
        % flux calculation    %
        %%%%%%%%%%%%%%%%%%%%%%%
        Fl_low  = B_left*uminmod;
        Fl_high = B_left*uminmod + 0.5*(1 - mu)*A_left*uminmod;

        Fr_low  = B_right*uminmod;
        Fr_high = B_right*uminmod + 0.5*(1 - mu)*A_right*uminmod;

        %%%%%%%%%%%%%%%%%%%%%%%
        % update u            %
        %%%%%%%%%%%%%%%%%%%%%%%
        C = minmod(r);
        F_left  = Fl_low + C.*(Fl_high - Fl_low);
        F_right = Fr_low + C.*(Fr_high - Fr_low);
        F = F_right - F_left;
        uminmod = uminmod - mu*F;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % superbee                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        % get ratios          %
        %%%%%%%%%%%%%%%%%%%%%%%
        r1 = A*usuperbee;
        r2 = B*usuperbee;
        r2(r2 == 0) = 1e-10;
        r = r1./r2;

        %%%%%%%%%%%%%%%%%%%%%%%
        % flux calculation    %
        %%%%%%%%%%%%%%%%%%%%%%%
        Fl_low  = B_left*usuperbee;
        Fl_high = B_left*usuperbee + 0.5*(1 - mu)*A_left*usuperbee;

        Fr_low  = B_right*usuperbee;
        Fr_high = B_right*usuperbee + 0.5*(1 - mu)*A_right*usuperbee;

        %%%%%%%%%%%%%%%%%%%%%%%
        % update u            %
        %%%%%%%%%%%%%%%%%%%%%%%
        C = superbee(r);
        F_left  = Fl_low + C.*(Fl_high - Fl_low);
        F_right = Fr_low + C.*(Fr_high - Fr_low);
        F = F_right - F_left;
        usuperbee = usuperbee - mu*F;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % van leer                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        % get ratios          %
        %%%%%%%%%%%%%%%%%%%%%%%
        r1 = A*uvanleer;
        r2 = B*uvanleer;
        r2(r2 == 0) = 1e-10;
        r = r1./r2;

        %%%%%%%%%%%%%%%%%%%%%%%
        % flux calculation    %
        %%%%%%%%%%%%%%%%%%%%%%%
        Fl_low  = B_left*uvanleer;
        Fl_high = B_left*uvanleer + 0.5*(1 - mu)*A_left*uvanleer;

        Fr_low  = B_right*uvanleer;
        Fr_high = B_right*uvanleer + 0.5*(1 - mu)*A_right*uvanleer;

        %%%%%%%%%%%%%%%%%%%%%%%
        % update u            %
        %%%%%%%%%%%%%%%%%%%%%%%
        C = vanleer(r);
        F_left  = Fl_low + C.*(Fl_high - Fl_low);
        F_right = Fr_low + C.*(Fr_high - Fr_low);
        F = F_right - F_left;
        uvanleer = uvanleer - mu*F;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % upwinding                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        % flux calculation    %
        %%%%%%%%%%%%%%%%%%%%%%%
        Fl_low  = B_left*uupwind;
        Fl_high = B_left*uupwind + 0.5*(1 - mu)*A_left*uupwind;

        Fr_low  = B_right*uupwind;
        Fr_high = B_right*uupwind + 0.5*(1 - mu)*A_right*uupwind;

        %%%%%%%%%%%%%%%%%%%%%%%
        % update u            %
        %%%%%%%%%%%%%%%%%%%%%%%
        C = zeros(size(r));
        F_left  = Fl_low + C.*(Fl_high - Fl_low);
        F_right = Fr_low + C.*(Fr_high - Fr_low);
        F = F_right - F_left;
        uupwind = uupwind - mu*F;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % exact solution                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:m
        if x(j) <= dt*t
            uexact(j) = u_l;
        else
            uexact(j) = u_r;
        endif
    endfor
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot solutions                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot(x, uminmod, 'r.-'); 
    plot(x, usuperbee, 'b.-'); 
    plot(x, uvanleer, 'g.-'); 
    plot(x, uupwind, 'c.-'); 
    plot(x, uexact, 'k.-'); 
    legend("minmod", "superbee", "van leer", "upwind", "exact");
    if t == 10
        print "-S640,480" prob4-1.png
    elseif t == timesteps/3
        print "-S640,480" prob4-2.png
    elseif t == timesteps/2
        print "-S640,480" prob4-3.png
    endif
    %savefig("fig%i.png", t);
    axis([a,b,u_r-0.1,u_l+0.1])
    drawnow
    clf;
endfor
input("Press any key to continue.");
