% Problem 5
% Marty Fuhry
% 1/20/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Measure the truncation error for finite difference approximations            %
% for the function                                                             %
%   u(x) = exp(sin(x))                                                         %
% over [-pi, pi].                                                              %
% Assume periodic boundary conditions.                                         %
%                                                                              %
% Left-Sided : u'(x)  = (u(x + dx) - u(x)) / dx                                %
% Right-Sided: u'(x)  = (u(x) - u(x - dx)) / dx                                %
% Centered   : u''(x) = (u(x + dx) - 2u(x) + u(x - dx)) / dx^2                 %
% Problem 4  : u'(x)  = (u(x - 2dx) - 6u(x - dx) + 3u(x) + 2u(x + dx)) / (6dx) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Octave is dumb; make sure it knows this isn't a function file
1;

% These will return the matrices for each finite difference method
% 1st-Order Left Sided
function [A] = fo_ls(m, dx)
    A = sparse(diag(-1*ones(m+2,1)) + diag(ones(m+1,1),1));
    A(m+2,2) = 1; % for BCs
    A *= 1/dx;
endfunction

% 1st-Order Right Sided
function [A] = fo_rs(m, dx)
    A = sparse(diag(1*ones(m+2,1))) + sparse(diag(-1*ones(m+1,1),-1));
    A(1,m+1) = -1; % for BCs
    A *= 1/dx;
endfunction

% 2nd-Order Centered
function [A] = so_c(m, dx)
    A = sparse(diag(ones(m+1,1),-1)) + sparse(diag(-2*ones(m+2,1))) + sparse(diag(ones(m+1,1), 1)); 
    A(1,m+1) = 1; % for BCs
    A(m+2,2) = 1; % for BCs
    A *= 1/(dx^2);
endfunction

% Problem 4
function [A] = p4(m, dx)
    A = sparse(diag(ones(m,1),-2));
    A += sparse(diag(-6*ones(m+1,1), -1));
    A += sparse(diag(3*ones(m+2,1)));
    A += sparse(diag(2*ones(m+1,1),1));
    A(1,m)   = 1;  % for BCs
    A(1,m+1) = -6; % for BCs
    A(2,m+1) = 1;  % for BCs
    A(m+2,2) = 2;  % for BCs
    A *= 1/(6*dx);
endfunction

clf; 
maxSize = 13; % If Octave is 64-bit and compiled with --enable-64, I could make this much bigger,
              % but there is a hard limit on the size of a vector in this version, so I can't make 
              % anything larger than 2^13. The machine doesn't even come close to running out of 
              % memory...

err_l = zeros(maxSize-3,1); % Error vectors
err_r = zeros(maxSize-3,1); 
err_c = zeros(maxSize-3,1);
err_4 = zeros(maxSize-3,1);

i = 1; % iterator for the errors

% Generate the approximation for each value of m
for m = 2.^[3:maxSize]
    dx = 2*pi / (m+1);      % size of delta x
    x = [-pi : dx : pi]';   % discretize the region

    u = exp(sin(x));                         % evaluate u(x)
    up  = cos(x).*exp(sin(x));               % u'(x)
    upp = (cos(x).^2 - sin(x)).*exp(sin(x)); % and u''(x)

    % 1st-Order Left-Sided
    A         = fo_ls(m, dx);
    up_l      = A * u;
    err_l(i)  = sqrt(dx)*norm(up_l - up, 2); % grid-point function 2-norm

    % 1st-Order Right-Sided
    A         = fo_rs(m, dx);
    up_r      = A * u;
    err_r(i)  = sqrt(dx)*norm(up_r - up, 2); % grid-point function 2-norm

    % 2nd-Order Centered
    A         = so_c(m, dx);
    up_c      = A * u;
    err_c(i)  = sqrt(dx)*norm(up_c - upp, 2); % grid-point function 2-norm

    % Problem 4
    A         = p4(m, dx);
    up_4      = A * u;
    err_4(i)  = sqrt(dx)*norm(up_4 - up, 2); % grid-point function 2-norm

    i += 1;

endfor

% Plot the errors
hold on;
title('Errors of Finite Difference Methods');
xlabel('m');
ylabel('Error (2-Norm)');
m = 2.^[3:maxSize];
loglog(m, err_l, 'r*-');
loglog(m, err_r, 'bo-');
loglog(m, err_c, 'g-');
loglog(m, err_4, 'k');
legend('Left-Sided', 'Right-Sided', 'Centered', 'Problem 4', 'location', 'southwest');
axis([2^3, 2^maxSize])

% This clearly shows the truncation error going to 0.
