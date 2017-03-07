function [post prob acti log_likelihood main mix errlog mix0 errlogKM] = em_gmm_v1(matrix,covartype,number_classes,varargin)
% em_gmm (deprec) Train a GMM model with EM algorithm and KMean initialisation.
%
% [POST PROB ACTI LLH LABELS MIX ERRLOG MIXKMEAN ERRLOGKMEAN] = em_gmm_v1(DATA,COVARTYPE,K,[PAR,VAL]) 
%
% Inputs:
%	DATA double array [Nd,Np]: the dataset, a collection of Np vectors of Nd dimensions.
% 	COVARTYPE: a string to set the type of Netlab GMM covariance type to use. It can be:
% 		- 'spherical'
% 		- 'diag'
% 		- 'full'
% 	K: an integer defining the number of components of the GMM
% 	PAR/VAL:
% 		- NiterKM
% 		- NrunsKM
% 		- NiterGM
% 		- NrunsGM
% 		- check_covar
% 		- debug
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2015-01-21 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

% This file is part of OBIDAM/PCM.
%     OBIDAM/PCM is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     OBIDAM/PCM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.


error('This function is no longer supported, please use em_gmm_v2')

%- Default parameters:
maxIterationsKMEAN = 100;
maxIterationsGMM   = 1000;
check_covar = 0;
debug = 0;

%- Load user options:
if nargin-3 > 0
	if mod(nargin-3,2) ~= 0
		error('Options must come in pairs: OPTION,VALUE')
	end% if 
	for in = 1 : 2 : nargin-3
		eval(sprintf('%s = varargin{in+1};',varargin{in}));
	end% for in	
	clear in
end% if

%- Variables
dimension = size(matrix,1);

%- Create a Gaussian mixture model structure
mix0 = gmm(dimension, number_classes, covartype);

%- Initialise GMM with a Kmean
options = foptions;
options(14) = maxIterationsKMEAN;	% Iterations of k-means in initialisation
[mix0, errlogKM] = gmminit(mix0, matrix', options);

%- Train GMM with EM

%-- Options
options = foptions;
options(1)  = 0;		% Prints out error values.
options(3)  = 0;		% No convergence level criteria
options(5)  = check_covar; % Ensure that covariances don't collapse
options(14) = maxIterationsGMM;		% Number of iterations.
options(19) = debug; % Debug by plot in live EM iterations

%-- Training
[mix, options, errlog] = gmmem(mix0, matrix', options, [], []);

%- Compute other usefull metrics for outputs
post = gmmpost(mix,matrix');
prob = gmmprob(mix,matrix');
acti = gmmactiv(mix,matrix');
log_likelihood = options(8);
[~,main] = max(post,[],2);


end %functionem_gmm