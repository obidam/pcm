function varargout = pcmprctilelatlon(GMMMASK,PCM,DPT,FIELD,M,varargin)
% pcmprctilelatlon Compute the percentile at each depth levels from a PCM and a given collection of gridded profiles
%
% Y = pcmprctilelatlon(MASK,PCM,DPT,DATA,M,[PAR,VAL])
%
% REQUIRED INPUTS:
% 	MASK(Ny,Nx): 0/1 mask to define profiles to classify in the lat/lon grid
%	PCM (struct): A Profile Classification Model structure
% 	DPT(Nz,1): Depth axis of DATA (used to call pcmpredict and get the activation values)
% 	DATA(Nz,Ny,Nx): A gridded collection of profiles
% 	M: scalar or vector of percent values between 0 and 100
% 
% OPTIONAL INPUTS as [PAR,VAL]: List of optional parameters/values: 
% 	REDUCED: 0 (default) / 1: Indicate if DATA is reduced or not. If it's reduced
% 		real data will be recomposed using MODEL info.
% 	WEIGHT(Ny,Nx,K): Possibly provide the profile weights to use to compute histograms.
% 		If not provided, DATA is classified using pcmpredict and ACTIVATIONS are used.
% 	CLASS: A list of integers to specify the class to determine
%
% OUTPUTS:
% 	Y: Percentiles of the values in DATA and CLASS
% 
% RQ:
% 	This is a simple wrapper of the pcmprctile function to handle gridded data
%
% See also: pcmprctile, pcmpdfzlatlon
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


%- Load user parameters:
rnin = 5;
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

%- Move profiles from Z,Y,X to Z,P:
c = GMMMASK'; % Y,X to X,Y
d = permute(FIELD,[1 3 2]); % Z,Y,X to Z,X,Y
ikeep = find(c(:)==1);
X = d(:,ikeep); % Z, P for pfoiles over the mask=1
[~, Np] = size(X);

if exist('WEIGHT','var')
	d = permute(WEIGHT,[3 2 1]); % Y,X,K to K,X,Y
	W = d(:,ikeep); % Z, P for profiles over the mask=1	
	W = W'; % because pcmprctile works with K,Np
	%- Classification:
	Y = pcmprctile(PCM,DPT,X,M,varargin{:},'WEIGHT',W);
else
	%- Classification:
	Y = pcmprctile(PCM,DPT,X,M,varargin{:});
end% if


%- Outputs
varargout(1) = {Y};

end %functionpcmprctilelatlon