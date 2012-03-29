% Problem 4 (a)
% Marty Fuhry
% 3/28/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a finite element code to solve the one-dimensional boundary value      %
% problem for u(x) on 0 < x < 1                                                %
%  u''(x) = e^(-10x), u(0) = 0, u(1) = 0.                                      %
% Use tent functions on 10 uniformly spaced elements.                          %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build the grid
h = 0.1;
x = [0:h:1]';
m = size(x)(1)-2;

K = zeros(m+2,m+2);
L = zeros(m+2,1);

% define the basis function over the canonical element
function [feval] = f(x)
    feval = exp(-10.*x);
endfunction

% build the element stiffness matrix
for j = 2:m+2
    K(j-1: j, j-1: j) += [1 -1; -1 1];
endfor

% boundary conditions
K(m+2,m+1) = 0;
K(1,2) = 0;

% build the load vector
F = [f(x)];
for j = 3:m+1
    L(j-1:j) += h/6. * [2*F(j-1) + F(j); F(j-1) + 2*F(j)];
endfor

% solve the steady state system
C = K \ L;

% and plot
plot(x,C);
drawnow
input("Press any key to continue.")
