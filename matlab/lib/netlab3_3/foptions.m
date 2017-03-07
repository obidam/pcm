function opt_vect = foptions()
% FOPTIONS Sets default parameters for optimisation routines
% For compatibility with MATLAB's foptions()
%
% options(1) = -1/0/1 : Print out error values
% options(14) : Nb of iterations
% options(7) = 1;    % Set width factor of RBF

% options(5) = 1;		% Use persistence
% options(7) = 50;	% Number of steps in trajectory.
% options(14) = nsamples;	% Number of Monte Carlo samples returned. 
% options(15) = 30;	% Number of samples omitted at start of chain.

% Copyright (c) Dharmesh Maniyar, Ian T. Nabney (2004)
  
opt_vect      = zeros(1, 19);
opt_vect(2:3) = 1e-4;
opt_vect(4)   = 1e-6;
opt_vect(16)  = 1e-8;
opt_vect(17)  = 0.1;