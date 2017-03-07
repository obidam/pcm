function DIST = pcmdist(MODEL,DPT,X)
% pcmdist Compute the Mahalanobis distance of a collection of profiles to PCM Gaussian Modes
%
% [] = pcmdist(PCM,DPT,DATA) HELP_TEXT_1ST_FORM
%
% REQUIRED INPUTS:
%
% OUTPUTS:
%
% EG:
%
% See Also: 
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-06-07 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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

%- Default parameters:
usemodelnorm = true;

%- Load user parameters:
if nargin > 3
    if mod(nargin-3,2) ~= 0
        error('Parameters must come in pairs: PAR,VAL');
    end% if
    for in = 1 : 2 : nargin-3
        eval(sprintf('%s = varargin{in+1};',varargin{in}));
    end% for in
    clear in
elseif nargin == 3
	% OK
else
	error('Bad number of parameters');
end% if

%- Interpolation (X->Xi): 
% Check if we need to interpolate the collection to be classified onto the PCM vertical grid
DPTmodel = MODEL.DPTmodel;

% Compare with the depth axis of the profiles to be classified (DPT) with the PCM depth axis:
[Nz, Np] = size(X);
if Nz == length(DPTmodel) & prod(DPT == DPTmodel) == 1
	% Similar axis, dont need to interpolate
	doINTERPz = false;
else
	% Different axis, we need to interpolate the input field
	doINTERPz = true;
end% if 

switch doINTERPz
	case 0 % No vertical interpolation
		Xi = X;
	case 1 % Run vertical interpolation
		% Possibly Create a mixed layer for the interpolation to work smoothly at the surface
		if (DPT(1)<0 & DPTmodel(1) == 0) & ~isnan(DPTmodel(1))
			DPT = [0; DPT];	
			X = cat(1,zeros(1,Np).*NaN,X);
			for ip = 1 : Np
				iz = find(~isnan(X(:,ip)),1,'first');
				X(1,ip) = X(iz,ip);
			end% for ip
		end% if
		% Interpolation
		if Np>1
			Xi = interp2(DPT,1:Np,X',DPTmodel,1:Np)';
		else
			Xi = interp1(DPT,X',DPTmodel)';			
		end% if 
end% switch 

fr = length(find(isnan(Xi(:))==1));
if fr > 1
	error(sprintf('I found NaNs into the working dataset (%i/%i),%s\n%s',...
		fr,numel(Xi(:)),...
		'I cannot classify these data.',...
		'This may be due to the model depth axis DPTmodel to be out of range of the DATA depth axis DPT'));
end% if 

%- Normalization (Xi->Xn):
switch MODEL.normalization
	case 0 % No norm
		Xn = Xi; 
	case 1 % Center/standardize
		switch usemodelnorm
			case 0 % Use DATA
				X_ave = nanmean(Xi,2);
				X_std = nanstd(Xi,[],2);
				Xn = (Xi - repmat(X_ave,[1 Np]))./repmat(X_std,[1 Np]);
			case 1 % Use model DATA
				Xn = (Xi - repmat(MODEL.X_ave,[1 Np]))./repmat(MODEL.X_std,[1 Np]);
		end% switch
	case 2 % Center only
		switch usemodelnorm
			case 0 % Use DATA
				X_ave = nanmean(Xi,2);
				Xn = (Xi - repmat(X_ave,[1 Np]));
			case 1 % Use model DATA
				Xn = (Xi - repmat(MODEL.X_ave,[1 Np]));
		end% switch
end% switch

%- Reduction (Xn->Xr):
% Reduce dimensions using GMM model PCA reduced dimensional basis
switch MODEL.doREDUCE
	case 0 % No reduction
		Xr = Xn;
	case 1 % Recution using PCA (from MODEL new space)
		Xr = MODEL.EOFs'*(Xn-repmat(MODEL.X_ref,[1 Np]));
end% switch 

%- Compute the distance:
% ACTI = gmmactiv(MODEL.mix,Xr');

switch MODEL.covarTYPE
	case 'full'
		for k = 1 : MODEL.K
			if 0
				diffs = Xr' - (ones(Np, 1) * MODEL.mix.centres(k, :));
				% We use the Cholesky decomposition of the covariance
				% matrix to speed up the computation:
				c = chol(MODEL.mix.covars(:, :, k));
				temp = diffs/c;
				DIST(:,k) = sqrt(sum(temp.*temp, 2));
			else
				% This is equivalent to with X [n * d]: 
				% X = x - repmat(mix.centres(j,:),[ndata 1]); % n * d
				% X = X'; % d * n, the correct dimension for Mahalanobis distance
				% z = mix.covars(:,:,j); % The squared d*d covariance matrix           
				% a(:,j) = diag(exp(-0.5*(X'*inv(z)*X))./normal./sqrt(det(z)))
				X = Xr - repmat(MODEL.mix.centres(k,:)',[1 Np]); % d * n
				z = MODEL.mix.covars(:,:,k); % The squared d*d covariance matrix 
				DIST(:,k) = diag(X'*inv(z)*X);
			end% if 
		end% for j
	otherwise
		error('Not implemented yet');
end% switch 


end %functionpcmdist