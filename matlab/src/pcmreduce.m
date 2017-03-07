function varargout = pcmreduce(MODEL,DPT,X,varargin)
% pcmreduce Reduce a data set using a given PCM
%
% [DATAr DATAn DATAi] = pcmreduce(PCM,DPT,DATA)
% 
% Reduce the data set DATA using the method and information from PCM
%
% REQUIRED INPUTS:
% 	PCM: A Profile Classification Model as output from PCMTRAIN
% 	DPT(Nz): The vertical axis of DATA
% 	DATA(Nz,Np): The dataset to reduce
%
% OPTIONAL INPUTS as [PAR,VAL]: List of optional parameters/values:
% 	'usemodelnorm': Determine which data to use during the normalization step
% 		true,1 (default):  Use the PCM stored data
% 		false,0: or use DATA to apply the normalization method
%
% See Also: pcmtrain, reduce_dimensions
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-05-20 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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
rnin = 3; % Number of required inputs
if nargin > rnin
    if mod(nargin-rnin,2) ~= 0
        error('Parameters must come in pairs: PAR,VAL');
    end% if
    for in = 1 : 2 : nargin-rnin
        eval(sprintf('%s = varargin{in+1};',varargin{in}));
    end% for in
    clear in
elseif nargin == rnin
	% OK
else
	error('Bad number of parameters');
end% if
clear rnin

%- Interpolation (X->Xi): 
% Check if we need to interpolate the collection to be reduced onto the PCM vertical grid
DPTmodel = MODEL.DPTmodel;

% Compare with the depth axis of the profiles to be reduced (DPT) with the PCM depth axis:
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
		Xi = interp2(DPT,1:Np,X',DPTmodel,1:Np)';
end% switch 

fr = length(find(isnan(Xi(:))==1));
if fr > 1
	error(sprintf('I found NaNs into the working dataset (%i/%i),%s\n%s',...
		fr,numel(Xi(:)),...
		'I cannot reduce these data.',...
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

%- Output
varargout(1) = {Xr};
varargout(2) = {Xn};
varargout(3) = {Xi};

end %functionpcmreduce