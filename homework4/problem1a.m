% Problem 1 (a)
% Marty Fuhry
% 3/27/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let u(x) = 1/(1 + 16x^2) over the interval [-1,1].                           %
% Use polynomial interpolation to approximate the derivative of u(x) over      %
%  x_j = -1 + j * dx                                                           %
% for dx = 2/m and m = 16.                                                     %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 16;
dx = 2/m;

x = [-1: dx :1];

u = 1./(1. + 16.*x.^2);
uprime = -(1. + 16.*x.^2).^(-2).*32.*x;

p = polyfit(x, u, m);
y = polyval(p, x);

pprime = polyderiv(p);
yprime = polyval(pprime, x);

hold on
plot(x, uprime)
plot(x, yprime)
drawnow

input("Press any key to continue.");
