function varargout = pcmpdfzlatlon(GMMMASK,PCM,DPT,FIELD,M,varargin)
% pcmpdfzlatlon Compute the PDF at each depth levels (PDFz) from a PCM and a given collection of gridded profiles
%
% [PDFz,X,PDFzO,PDFzM] = pcmpdfzlatlon(MASK,PCM,DPT,DATA,M,[PAR,VAL])
%
% REQUIRED INPUTS:
% 	MASK(Ny,Nx): 0/1 mask to define profiles to classify in the lat/lon grid
%	MODEL (struct): A Profile Classification Model structure
% 	DPT(Nz,1): Depth axis of DATA (used to call pcmpredict and get the activation values)
% 	DATA(Nz,Ny,Nx): A gridded collection of profiles
% 	M: if M is a scalar, uses M bins, but if M is a vector, returns 
% 		the PDFz of DATA among bins with centers specified by M.
% 
% OPTIONAL INPUTS as [PAR,VAL]: List of optional parameters/values: 
% 	REDUCED: 0 (default) / 1: Indicate if DATA is reduced or not. If it's reduced
% 		real data will be recomposed using MODEL info.
% 	WEIGHT(Np,K): Possibly provide the profile weights to use to compute histograms.
% 		If not provided, DATA is classified using pcmpredict and ACTIVATIONS are used.
% 	CLASS: A list of integers to specify the class to determine
%
% OUTPUTS:
% 	PDFz(Nz,X,K): PDFz for each class of the model
% 	X: Also returns the position of the bin centers
% 	PDFzO(Nz,M): Observed PDFz from the entire collection
% 	PDFzM(Nz,M): Model PDFz from the prior weighted sum of PDFz
%
% EG:
%
% 
% See Also: pcmpdfz, pcmprctilelatlon
% 
% Reference:
% 	Maze et al, Prog. Oc. (2016): "Coherent heat structures revealed by unsupervised classification 
% 	of Argo temperature profiles in the North Atlantic Ocean".
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-04-15 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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
[PDFz,X,PDFzO,PDFzM] = pcmpdfz(PCM,DPT,X,M,varargin{:});

%- Outputs
varargout(1) = {PDFz};
varargout(2) = {X};
varargout(3) = {PDFzO};
varargout(4) = {PDFzM};



end %functionpcmpdfzlatlon