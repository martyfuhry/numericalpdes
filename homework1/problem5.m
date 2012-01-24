% Problem 5
% Marty Fuhry
% 1/20/2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure the truncation error for finite difference approximations         %
% for the function                                                          %
%   u(x) = exp(sin(x))                                                      %
% over [-pi, pi].                                                           %
% Assume periodic boundary conditions.                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u'(x) = (u(x + dx) - u(x)) / dx
% u'(x) = (u(x) - u(x - dx)) / dx
% u''(x) = (u(x + dx) - 2u(x) + u(x - dx)) / dx^2
% u'(x) = (u(x - 2dx) - 6u(x - dx) + 3u(x) + 2u(x + dx)) / dx

clf;

% 1st-Order Left Sided
function [A] = fo_ls(m, dx)
    A = sparse(diag(-1*ones(m+2,1)) + diag(ones(m+1,1),1));
    A(m+2,2) = 1;
    A *= 1/dx;
endfunction

% 1st-Order Right Sided
function [A] = fo_rs(m, dx)
    A = sparse(diag(1*ones(m+2,1)) + diag(-1*ones(m+1,1),-1));
    A(1,m+1) = -1;
    A *= 1/dx;
endfunction

% 2nd-Order Centered
function [A] = so_c(m, dx)
    A = sparse(diag(ones(m+1,1),-1) + diag(-2*ones(m+2,1)) + diag(ones(m+1,1), 1)); 
    A(1,m+1) = 1;   % boundary conditions for u_0
    A(m+2,2) = 1;   % boundary conditions for u_m+1
    A *= 1/(dx^2);
endfunction

% Problem 4
function [A] = p4(m, dx)
    A = sparse(diag(ones(m,1),-2));
    A += diag(-6*ones(m+1,1), -1);
    A += diag(3*ones(m+2,1));
    A += diag(2*ones(m+1,1),1);
    A(1,m) = 1;
    A(1,m+1) = -6;
    A(2,m+1) = 1;
    A(m+2,2) = 2;
    A *= 1/(6*dx);
endfunction

for m = [80]
    dx = 2*pi / (m+1);      % size of delta x
    x = [-pi : dx : pi]';   % discretize the region
    u = exp(sin(x));        % evaluate u

    % 1st-Order Left-Sided
    %A = fo_ls(m, dx);
    %up_l = A * u;

    % 1st-Order Right-Sided
    %A = fo_rs(m, dx);
    %up_r = A * u;

    % 2nd-Order Centered
    %A = so_c(m, dx);
    %upp_c = A * u;

    % Problem 4
    A = p4(m, dx);
    up_4 = A * u;


    % Plot the errors
    up  = cos(x).*exp(sin(x));
    upp = (cos(x).^2 - sin(x)).*exp(sin(x));

    hold on;

    %plot(x,upp_c,'r');
    %plot(x,upp,'b');

    plot(x,up_4,'r');
    plot(x,up,'b');

endfor
