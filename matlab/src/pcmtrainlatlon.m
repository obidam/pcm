function varargout = pcmtrainlatlon(GMMMASK,FIELD,K,COVARTYPE,DPTmodel,varargin)
% pcmtrainlatlon Train a Profile Classification Model (PCM) based on GMM onto a depth/lat/lon gridded field
% 
% [PCM DATAi DATAn DATAr] = pcmtrainlatlon(MASK,DATA,K,COVARTYPE,DPTmodel,[PAR,VAL])
%
% REQUIRED INPUTS:
% 	MASK(Ny,Nx): 0/1 mask to define profiles to classify on the lat/lon grid
% 	DATA(Nz,Ny,Nx): The 3D field to classify with a GMM, Ny.Nx profiles with Nz depth levels.
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
% 	maxvar: The maximum variance to retained during compression with PCA if doREDUCE==1
% 	
% OUTPUTS:
% 	PCM: a Profile Classification Model (PCM) structure
% 	DATAi: DATA after the interpolation step (if interp was not necessary, this is DATA)
% 	DATAn: DATA after normalization step (if normalization=0, this is DATAi)
% 	DATAr: DATA after the reduction step (if doREDUCE=0, this is DATAn)
% 
% Workflow of the function:
% 	Move profiles from Z,Y,X to Z,P (according to MASK values)
% 	Interpolate profiles if necessary (according to DPT and DPTmodel)
% 	Normalize profiles if necessary (according to normalization value)
% 	Reduce dimensions if necessary (if doREDUCE true and according to maxvar value)
% 	Train a GMM (according to K and COVARTYPE values)
% 	Move normalized and reduced fields back to the lat/lon grid
% 	Outputs (format the PCM structure)
%
% Rq:
% 	- DPTmodel and DPT axis must be negative, oriented from the surface toward the bottom
% 	- This is a simple wrapper of the pcmtrain function to handle gridded data
% 
% See also: pcmpredictlatlon, pcmtrain
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2015-11-06 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)
% Revised: 2016-04-15 (G. Maze) Delegate the training to the pcmtrain function

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


%- Load user parameters:
if nargin > 5
    if mod(nargin-5,2) ~= 0
        error('Optional parameters must come in pairs: PAR,VAL');
    end% if
    for in = 1 : 2 : nargin-5
        eval(sprintf('%s = varargin{in+1};',varargin{in}));
    end% for in
    clear in
elseif nargin == 5
	% OK
else
	error('Bad number of parameters');
end% if

%- Move profiles from DPT,LAT,LON to DPT,PROF:
[Nz, Ny, Nx] = size(FIELD);
c = GMMMASK'; % X, Y
d = permute(FIELD,[1 3 2]); % Z, X, Y
ikeep = find(c(:)==1);
X = d(:,ikeep);
[~, Np] = size(X);

%- Train model:
[MODEL, Xi, Xn, Xr] = pcmtrain(X,K,COVARTYPE,DPTmodel,varargin{:});

%- Move back onto the lat/lon grid:

% Init
Fi = zeros(size(Xi,1),Nx,Ny)*NaN;
Fn = zeros(size(Xn,1),Nx,Ny)*NaN;
Fr = zeros(size(Xr,1),Nx,Ny)*NaN;

% Fill
Fi(:,ikeep) = Xi;
Fn(:,ikeep) = Xn; 
Fr(:,ikeep) = Xr; 

% Format
Fi = permute(Fi,[1 3 2]); % Nz,Ny,Nx
Fn = permute(Fn,[1 3 2]); % Nz,Ny,Nx
Fr = permute(Fr,[1 3 2]); % Ny,Nx,Neofs

%- Outputs
% Un-used fields are simply not added to the model structure
MODEL.readme = 'This PCM was created using pcmtrainlatlon.m';

varargout(1) = {MODEL};
varargout(2) = {Fi};
varargout(3) = {Fn};
varargout(4) = {Fr};


end %functionpcmtrainlatlon