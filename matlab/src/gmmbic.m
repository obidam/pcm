function [Y Y1 Y2 n] = gmmbic(mix,nllh,x,varargin)
% gmmbic Return the BIC from a GMM and a dataset
%
% BIC = gmmbic(MIX,NLLH,X) Return the BIC from a GMM and a dataset
%
% 	BIC = 2*NLLH  +  NumParams*Log(NumObs)
% 
% [Nd NumObs] = size(X);
%
% PCM is Profile Classification Modelling
% Copyright (C) 2015-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2015-12-07 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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


%- Get nb of observations and dimensions:
[Nd Nobs] = size(x);
K = mix.ncentres;

if Nd ~= mix.nin
	error('The GMM and dataset number of dimensions do not match !')
end% if 

%- Compute the number of free parameters to estimate:
% Nb of priors:
n = K-1;

% Add Nb of Gaussian centers:
n = n + K*Nd;

% Add Nb of Gaussian covariance matrix elements:
switch mix.covar_type
	case 'full'
		n = n + K*Nd*(Nd-1)/2;
	case 'diag'
		n = n + K*Nd;
	case 'spherical'
		n = n + K;
	otherwise
		error('Unknown model type of covariance matrix')		
end% switch 

%- Compute BIC = -2*LLH  +  NumParams*Log(NumObs)
%Y  = 2*nllh + n*log(Nobs);
Y1 = 2*nllh;
Y2 = n*log(Nobs);
Y = Y1 + Y2;

end %functiongmmbic