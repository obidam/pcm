function [POST, LABELS, ACTI, PROB] = pcmpredictlatlon(GMMMASK,MODEL,DPT,FIELD,varargin)
% pcmpredictlatlon Classify a lat/lon gridded field with Profile Classification Model based on GMM 
%
% [POST LABELS ACTI PROB] = pcmpredictlatlon(MASK,PCM,DPT,FIELD) 
%
% Inputs:
% 	MASK(Ny,Nx): 0/1 mask to define profiles to classify in the lat/lon grid
% 	PCM: Profile Classification Model (PCM) structure:
% 		- mix: The netlab GMM model structure
% 		- DPT: The vertical axis GMM has been trained on
% 		- EOFs: The eigenvectors used to compress data
% 		- Xref: The eigenvectors reference
% 	DPT(Nz,1): The vertical axis of FIELD (required only if the PCM uses vertical interpolation)
% 	FIELD(Nz,Ny,Nx): The 3D field to apply GMM on
% PAR,VAL: List of optional parameters/values: 
% 	'usemodelnorm':
% 		true,1:  Will use the model data for the normalization
% 		false,0: Will use FIELD data to apply the normalization method
%
% Rq:
% 	FIELD data are normalized along the sample axis.
% 	DPT axis is considered negative from surface to bottom.
%
% See also: pcmtrainlatlon
% 
% PCM is Profile Classification Modelling
% Copyright (C) 2015-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2015-10-31 ((G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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


%- Move profiles from Z,Y,X to Z,P:
c = GMMMASK'; % Y,X to X,Y
d = permute(FIELD,[1 3 2]); % Z,Y,X to Z,X,Y
ikeep = find(c(:)==1);
X = d(:,ikeep); % Z, P for pfoiles over the mask=1
[~, Np] = size(X);

%- Classification:
[post labels acti prob] = pcmpredict(MODEL,DPT,X,varargin{:});

%- Move results back on the original Y,X grid:
% Init
[Nz Ny Nx] = size(FIELD);
K = MODEL.mix.ncentres;	
POST = zeros(K,Nx,Ny)*NaN;
LABELS = zeros(Nx,Ny)*NaN;
ACTI = zeros(K,Nx,Ny)*NaN;
PROB = zeros(Nx,Ny)*NaN;
% Fill
POST(:,ikeep) = post'; 
LABELS(ikeep) = labels';
ACTI(:,ikeep) = acti'; 
PROB(ikeep) = prob';
% Format
POST = permute(POST,[3 2 1]); % Ny,Nx,K
LABELS = LABELS'; % Ny,Nx
ACTI = permute(ACTI,[3 2 1]); % Ny,Nx,K
PROB = PROB'; % Ny,Nx

end %functionpcmpredictlatlon