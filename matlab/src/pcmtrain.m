function varargout = pcmtrain(X,K,COVARTYPE,DPTmodel,varargin)
% pcmtrain Train a Profile Classification Model (PCM) based on GMM onto a collection of profiles
%
% [PCM DATAi DATAn DATAr] = pcmtrain(DATA,K,COVARTYPE,DPTmodel,[PAR,VAL])
%
% REQUIRED INPUTS:
% 	DATA(Nz,Np): The field to classify with a GMM, Np profiles with Nz depth levels.
% 	K: An integer to indicate the number of class to use.
% 	COVARTYPE: A string defining the shape of the Gaussian class covariance matrix:
% 		- 'spherical': Identity, a single standard deviation is used
% 		- 'diag': Diagonal matrix, anisotrop / orthogonal
% 		- 'full': Plain matrix, anisotrop / non-orthogonal
% 	DPTmodel(Nz,1): The vertical axis the PCM will work on. If the optional argument DPT is not
% 		provided, DPTmodel is not used and is simply added to the PCM structure.
% 		Rq: It will be mandatory in the prediction step (using pcmpredict or pcmpredictlatlon)
% 		or to save the model on disk (pcmsave).
% 
% OPTIONAL INPUTS as [PAR,VAL]: List of optional parameters/values: 
% 	DPT(Nz,1): The vertical axis of DATA. If it is different from DPTmodel, profiles will be
% 		interpolated onto DPTmodel using a simple 2D linear interpolation.
% 	normalization: An integer to choose the normalization method: 
% 		0: Nothing is done
% 		1 (default): Centered and standardized at each depth levels
% 		2: Centered at each depth levels
% 	doREDUCE: An integer to indicate if the vertical dimension should be reduced or not
% 		0: No reduction
% 		1 (default): Reduction is performed using Principal Components Method (PCA)
% 	maxvar: The maximum variance (in %) to be retained during compression with PCA if doREDUCE==1
% 	
% OUTPUTS:
% 	PCM: a Profile Classification Model (PCM) structure
% 	DATAi: DATA after the interpolation step (if interp was not necessary, this is DATA)
% 	DATAn: DATA after normalization step (if normalization=0, this is DATAi)
% 	DATAr: DATA after the reduction step (if doREDUCE=0, this is DATAn)
% 
% Workflow of the function:
% 	Interpolate profiles if necessary (according to DPT and DPTmodel)
% 	Normalize profiles if necessary (according to normalization value)
% 	Reduce dimensions if necessary (if doREDUCE true and according to maxvar value)
% 	Train a GMM (according to K and COVARTYPE values)
% 	Outputs (format the PCM structure)
%
% Rq:
% 	DPTmodel and DPT axis must be negative, oriented from the surface toward the bottom.
%
% See also: pcmpredict, pcmtrainlatlon
% 
% Reference:
% 	Maze et al, Prog. Oc. (2016): "Coherent heat structures revealed by unsupervised classification 
% 	of Argo temperature profiles in the North Atlantic Ocean".
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-04-14 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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
doINTERPz = false;
doREDUCE  = true;
maxvar    = 99.9;
normalization = 1;

%- Load user parameters:
if nargin > 4
    if mod(nargin-4,2) ~= 0
        error('Optional parameters must come in pairs: PAR,VAL');
    end% if
    for in = 1 : 2 : nargin-4
        eval(sprintf('%s = varargin{in+1};',varargin{in}));
    end% for in
    clear in
elseif nargin == 4
	% OK
else
	error('Bad number of parameters');
end% if

%- Size of the training set X:
[Nz, Np] = size(X);

%- Interpolation (X->Xi):
% Check if we need to interpolate the training set:
if exist('DPT','var')
	if Nz == length(DPTmodel) & prod(DPT == DPTmodel) == 1
		% Similar axis, dont need to interpolate
		doINTERPz = false;
	else
		% Different axis, we need to interpolate the inpute field
		doINTERPz = true;
	end% if
end% if 

switch doINTERPz
	case 0 %-- No interpolation
		Xi = X;
		if Nz ~= length(DPTmodel)
			error(sprintf('DATA has %i depth levels compared to %i in DPTmodel !',Nz,length(DPTmodel)))
		end% if
	case 1 %-- Linear interpolation
		DPT = DPT(:);
		DPTmodel = DPTmodel(:);
		if Nz ~= length(DPT)
			error(sprintf('DATA has %i depth levels compared to %i in DPT !',Nz,length(DPT)))
		end% if
		% Possibly Create a mixed layer for the interpolation to work smoothly at the surface		
		if DPT(1)<0 & DPTmodel(1) == 0
			DPT = [0; DPT];	
			X = cat(1,zeros(1,Np).*NaN,X);
			for ip = 1 : Np
				iz = find(~isnan(X(:,ip)),1,'first');
				X(1,ip) = X(iz,ip);
			end% for ip
		end% if
		% Interp this:
		Xi = interp2(DPT,1:Np,X',DPTmodel,1:Np)';
end% switch 

%- Normalization (Xi->Xn):
switch normalization
	case 0 %-- No norm
		Xn = Xi;
	case 1 %-- Center/standardize
		X_ave = nanmean(Xi,2);
		X_std = nanstd(Xi,[],2);
		Xn = (Xi - repmat(X_ave,[1 Np]))./repmat(X_std,[1 Np]);
	case 2 %-- Center only
		X_ave = nanmean(Xi,2);
		Xn = (Xi - repmat(X_ave,[1 Np]));
end% switch 

%stophere

%- Reduction (Xn->Xr):
switch doREDUCE
	case 0 	%-- No reduction
		Xr = Xn;
	case 1 	%-- Reduce dimension using PCA
		[Xr,EOFs,V,L,X_ref,Xp,Xrn,EOFsunit] = reduce_dimensions(Xn,maxvar); 
		% Xr: X reduced, Xp: reduced X re-projected on original data space
		maxvar = sum(V)*100; % Effective variance retained by the compression in %	
end% switch 

%- Training: 
% Train a GMM with Netlab throughout our customized workflow
[~,~,~,llh,~,mix] = em_gmm_v2(Xr,K,'covartype',COVARTYPE);

%- Outputs
% Un-used fields are simply not added to the model structure

MODEL = struct();
MODEL.Np = Np; % Size of the training set
MODEL.readme = 'This PCM was created using pcmtrain.m';

% Interpolation:
MODEL.DPTmodel  = DPTmodel;

% Normalisation:
MODEL.normalization = normalization;
switch normalization
	case 0
	case 1
		MODEL.X_ave = X_ave;
		MODEL.X_std = X_std;
	case 2
		MODEL.X_ave = X_ave;
end% switch 

% Reduction:
MODEL.doREDUCE = doREDUCE;
if doREDUCE
	MODEL.maxvar = maxvar;
	MODEL.EOFs = EOFs;
	MODEL.X_ref = X_ref;
	MODEL.V = V; % V [1,Nc] : The fraction of variance explained by individual axis of the new space.
%	MODEL.L = L; % L [1,Nc] : Principal component variances, ie eigenvalues of the covariance matrix.
end% if 

% Trained GMM:
MODEL.K = K;
MODEL.covarTYPE = COVARTYPE;
MODEL.mix = mix;
%MODEL.llh = llh; % This is wrong because em_gmm_v2 return in fact the negative llh
MODEL.llh = -llh; % Fix sign issue
MODEL.score = -llh/Np; % This is the sample-mean llh, compatible with scikit-learn "score" per-sample average 
                       % log-likelihood of the given data X
MODEL = orderfields(MODEL);

varargout(1) = {MODEL};
varargout(2) = {Xi};
varargout(3) = {Xn};
varargout(4) = {Xr};

end %functionpcmtrain