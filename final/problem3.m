% Problem 3 
% Marty Fuhry
% 4/15/2012
% Compiled and ran using GNU Octave, version 3.2.4 configured for "x86_64-pc-linux-gnu".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design & implement a finite element method to solve                          %
%  u_t = (kappa u_x)_x, u(0) = 1, u(1) = 1                                     %
% where                                                                        %
%  kappa(x) = 0.1 + 0.5(1 + tanh(6x - 3)                                       %
% and                                                                          %
%  u(x,0) = sech(200x - 100).                                                  %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build the grid
m = 128;
a = 0;
b = 1;

h = (b - a) / (m - 1);

x = [a:(b - a) / (2*m + 2):b]';
k = 0.1 + 0.5*(1 + tanh(6*x - 3));

K = zeros(2*m+3,2*m+3);
M = zeros(2*m+3,2*m+3);

dt = 0.000002;

% build the element stiffness matrix 
%   K_ij = (phi_i, phi_j) + dt * (phi_i', phi_j')_kappa
for j = 3:2*m+3
    % first add (phi_i, phi_j)
    K(j-2: j, j-2: j) += h * [2/15,  1/15, -1/30;
                              1/15,  8/15, 1/15;
                             -1/30, 1/15, 2/15;]; 

    % interpolate kappa(x) = k_j-1 * N_j-1 + k_j-1/2 N_j-1/2 + k_j N_j

    % k_j-1 * N_0
    K(j-2: j, j-2: j) += -k(j-2) * 2 * dt / h * [37/60 , -11/15, 7/60 ;
                                                 -11/15, 4/5   , -1/15;
                                                 7/60  , -1/15 , -1/20;];

    % k_j-1/2 * N_1
    K(j-2: j, j-2: j) += -k(j-1) * 2 * dt / h * [3/5  , -8/15, -1/15;
                                                 -8/15, 16/15, -8/15;
                                                 -1/15, -8/15, 3/5  ;];

    % k_j * N_2
    K(j-2: j, j-2: j) += -k(j) * 2 * dt / h * [-1/20, -1/15 , 7/60;
                                               -1/15, 4/5   , -11/15;
                                               7/60 , -11/15, 37/60];
endfor

% build the other matrix 
%   M_i,j = (phi_i, phi_j)
for j = 3:2*m+3
    M(j-2: j, j-2: j) += h * [2/15, 1/15, -1/30;
                              1/15, 8/15, 1/15;
                             -1/30, 1/15, 2/15;]; 
endfor


% initial conditions
C = sech(200*x - 100);

% boundary conditions 
K(1,1) = 1;
K(1,2) = 0;
K(1,3) = 0;
K(2*m+3,2*m+3) = 1;
K(2*m+3,2*m+2) = 0;
K(2*m+3,2*m+1) = 0;
M(1,1) = 1;
M(1,2) = 0;
M(1,3) = 0;
M(2*m+3,2*m+3) = 1;
M(2*m+3,2*m+2) = 0;
M(2*m+3,2*m+1) = 0;

C(1) = 1;
C(2*m+3) = 1;

for t = 0:0.001/dt
    % solve the system M * C_n+1 = K * C_n
    C = M \ (K * C);

    % and plot
    plot(x,C);

    legend("Finite Element Approximation");
    title("Problem 3");

    axis([a,b,0,1]);
    drawnow
endfor
%print("problem3.png")
input("Press any key to continue.");
