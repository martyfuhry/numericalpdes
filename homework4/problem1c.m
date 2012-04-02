% Problem 1 (c)
% Marty Fuhry
% 3/27/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot both interpolating polynomials on a fine uniformly spaced grid with     %
% m = 128 points. Discuss what you see (see paper).                            %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 16;
dx = 2/m;

x1 = [-1: dx :1];
x2 = cos([0:m]*pi / m);

x = [-1: 2/128 : 1];
u = 1./(1. + 16.*x.^2);
u1 = 1./(1. + 16.*x1.^2);
u2 = 1./(1. + 16.*x2.^2);

p1 = polyfit(x1, u1, m);
y1 = polyval(p1, x);

p2 = polyfit(x2, u2, m);
y2 = polyval(p2, x);

hold on
plot(x, u, 'r-.')
plot(x, y1, 'b-*')
plot(x, y2, 'k-x')

legend("Exact", "Uniform", "cos(j*pi/m)");
print("problem1b.png")
drawnow

input("Press any key to continue.");
