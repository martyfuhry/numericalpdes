% Problem 4 (b)
% Marty Fuhry
% 3/28/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a finite element code to solve the one-dimensional boundary value      %
% problem for u(x) on 0 < x < 1                                                %
%  u''(x) = e^(-10x), u(0) = 0, u(1) = 0.                                      %
% Use tent functions on 10 uniformly spaced elements.                          %
% Measure the error in (a) using the 2-norm. Can you reduce the error (by      %
% how much?) by using 10 non-uniformly spaced elements?                        %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build the grid
h = 0.1;

% uniform x; Error = 0.001243
x = [0:h:1]';
m = size(x)(1)-2;

% manually refined grid; Error = 0.001600
%x = [0, .13, .18, .23, .26, .28, .35, .5, .7, 1]';

% chebyshev nodes: Error = 0.001136
x = 0.5 + 0.5*cos((2*[m+2:-1:1] - 1)/(2*(m+2)) * pi)';
x(1) = 0;
x(m+2) = 1;

% calculate the h_j's
h = zeros(m+1,1);
for j = 1:m+1
    h(j) = x(j+1) - x(j);
endfor

K = zeros(m+2,m+2);
L = zeros(m+2,1);

% exact solution
u = (exp(-10*x) + (1 - exp(-10))*x - 1)/100;

% define the basis function over the canonical element
function [feval] = f(x)
    feval = -exp(-10.*x);
endfunction

% build the element stiffness matrix
for j = 2:m+2
    K(j-1: j, j-1: j) += 1/h(j-1) * [1 -1; -1 1];
endfor

% build the load vector
F = [f(x)];
for j = 2:m+2
    L(j-1:j) += h(j-1)/6. * [2*F(j-1) + F(j); F(j-1) + 2*F(j)];
endfor

% for the boundary conditions
K(m+2,m+1) = 0;
K(m+2,m+2) = 1;
K(1,2) = 0;
K(1,1) = 1;
L(1) = 0;
L(m+2) = 0;

% solve the steady state system
C = K \ L;

% and plot
hold on
plot(x,C,"b-.");
plot(x,u,"r-*");

legend("Finite Element Approximation", "Exact");
title("Problem 4 (b)");
print("problem4b.png");

printf("Error = %f\n", norm(u - C,2));
drawnow;
input("Press any key to continue.");
