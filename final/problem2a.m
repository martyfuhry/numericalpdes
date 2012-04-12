% Problem 2 (a)
% Marty Fuhry
% 4/11/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
alpha = 1;
m = 128;
a = -20;
b = 20;
timesteps = 1000;

k1 = [0:m/2]';
k2 = [-m/2 + 1: -1]';
k = [k1; k2];

dx = (b - a) / (m - 1);    % space discretization
dt = dx/(m^4/16); % not right either
x  = [a:dx:b]';      % grid

U = exp(-x.^2);
A = -diag(ones(m-1,1),-1) + diag(ones(m-1,1),1);
A(1,m) = -1;
A(m,1) = 1;

A *= 1/(4 * dx);

for t = 1:timesteps        % loop through the timesteps
    % RK4 time differencing

    % this is constant speed advection
    %k1 = dt * ifft((-i*k + k.^2 - k.^4).*fft(U));
    %k2 = dt * ifft((-i*k + k.^2 - k.^4).*fft(U + 0.5*k1));
    %k3 = dt * ifft((-i*k + k.^2 - k.^4).*fft(U + 0.5*k2));
    %k4 = dt * ifft((-i*k + k.^2 - k.^4).*fft(U + k3));

    % this works correctly
    k1 = dt* (-A * (U.^2))            + dt * ifft((k.^2 - k.^4).*fft(U));
    k2 = dt* (-A * ((U + 0.5k1).^2)) + dt * ifft((k.^2 - k.^4).*fft(U + 0.5*k1));
    k3 = dt* (-A * ((U + 0.5*k2).^2)) + dt * ifft((k.^2 - k.^4).*fft(U + 0.5*k2));
    k4 = dt* (-A * ((U + 1.0*k3).^2)) + dt * ifft((k.^2 - k.^4).*fft(U + 1.0*k3));

    % compute the result for this time step
    U += 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    plot(x,U);             % plot the result
    axis([a,b,-1,1]);
    drawnow;
endfor
