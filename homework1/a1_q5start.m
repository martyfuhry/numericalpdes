% Analyze truncation error in different finite difference approximations to du/dx
% Note, if this is your first time using Matlab, you should review
% the quickstart guide available on LEARN.  

gridsizes = 2.^(3:7); % Vector of various grid sizes

for m = gridsizes  % loop over the various m values and compute the
                   % truncation errors for each method

  dx = 
  x  = 
  u  = 
  uprime = 

  e  = ones(m,1);
  em = ones(m-1,1);
  emm= ones(m-2,1);

  % 1. 1st order right-sided difference:
  % Note: you will have better performance if you use the "sparse"
  % command to build your matrices -- type "help sparse" for syntax.
  D1 = diag(em,1)-diag(e,0); D1 = D1/dx;  % Diff. matrix
  error1 = norm(D1*u-uprime,inf);         % Truncation error

  loglog(m,[error1],'.','markersize',15.), hold on

end

legend('right-sided 1st order')
xlabel 'm'; ylabel 'error'
title  'Convergence of finite differences';