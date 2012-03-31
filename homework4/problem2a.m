% Problem 2 (a)
% Marty Fuhry
% 3/27/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the function                                                        %
%  u(x) = cos(k_max * x)                                                       %
% on the grid xj = j * dx for j = 1..m, dx = 2*pi/m, k_max = m/2.              %
% Compute by hand the approximation to u'(x_j) obtained by the pseudospectral  %
% method, and compare this to the exact value for u'(x_j). Why should we neg-  %
% lect the contributions                                                       %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 12;
dx = 2*pi/m;
kmax = m/2;

x = [dx:dx:2*pi]';

u = cos(kmax * x);

% for octave's stupid fourier coefficient indexing
k1 = [0:m/2]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

uprimeexact = kmax * (-sin(kmax * x));
uprimeps = ifft(i*k.*fft(u))


%ifft(i*k.*fft(u))
%hold on
%plot(x, uprimeexact, 'b-')
%plot(x,uprimeps, 'r*')

%drawnow

%input("Press any key to continue.");

