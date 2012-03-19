% Problem 1 (a)
% Marty Fuhry
% 3/15/2011
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the shallow water equations                                         %
%  u_t + n_x = 0                                                               %
%  n_t + u_x = 0                                                               %
% on the interval 0 < x < 1.                                                   %
% Use Lax Wendroff                                                             % 
%  u_j^n+1 = u_j^n - 1/(2dx) A (u_j+1^n - u_j-1^n)                             %
%            + dt^2/(dx^2) A (u_j+1^n - 2u_j^n + u_j-1^n)                      %
% with m=128 grid points and the intial condition                              %
%  n(x, 0) = sech(40(x - 0.5))                                                 %
%  u(x, 0) = 0                                                                 %
% with boundary condtions                                                      %
%  u(0, t) = 0                                                                 %
%  n(1, t) = 0                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up some PDE variables
m = 128;
a = 0;
b = 1;

% discretize the grid
dx = (b - a) / (m - 1);
x = [a-dx:dx:b+dx]';

% keep the courant number 1 so that there is no dispersion
dt = dx;
timesteps = 1/dt; % go until t = 1
mu = dt/dx;

% define the FD matrices
A = diag(ones(m+1,1),1) - diag(ones(m+1,1),-1);
B = diag(ones(m+1,1),1) - 2*diag(ones(m+2,1)) + diag(ones(m+1,1),-1);

% get the change of basis matrix P
M = [0 1; 1 0];
[P, L] = eig(M); % note: P is orthonormal because M is symmetric

% initial conditions
n = sech(40*(x - 0.5)); 
u = zeros(m+2,1);

% transform to the eigenbasis
W = P*[u'; n'];

w1 = W(1,:)';
w2 = W(2,:)';

% start the FD method
for t = 1:ceil(timesteps)
    % set left boundary conditions wrt the eigenbasis
    W = P*[0, n(2)]';
    w1(1) = W(1);
    w2(1) = W(2);

    % set right boundary conditions wrt the eigenbasis
    W = P*[u(m+1), 0]';
    w1(m+2) = W(1);
    w2(m+2) = W(2);

    % peform the FD method over the decoupled equations
    w1 = w1 + ( mu/2*A + mu^2/2*B)*w1;
    w2 = w2 + (-mu/2*A + mu^2/2*B)*w2;

    % transform back to the original basis
    W = [w1'; w2'];
    U = P'* W; 
    u = U(1,:)';
    n = U(2,:)';

    % plot u and n
    hold on
    plot(x,u, 'r'); 
    plot(x,n, 'b'); 
    legend("u", "n")
    axis([a,b,-1,1])
    drawnow;
    if t == 10
        print "-S640,480" prob1-1.png
    elseif t == ceil(timesteps/3)
        print "-S640,480" prob1-2.png
    elseif t == ceil(timesteps/2)
        print "-S640,480" prob1-3.png
    elseif t == ceil(timesteps)
        print "-S640,480" prob1-4.png
    endif
    clf;
endfor
input("Press any key to continue.");
