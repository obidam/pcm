function Xn = pcollecnorm(ARG1,Xi,varargin)
% pcollecnorm Normalize a collection of profiles
%
% Usage 1:
% 
% DATAn = pcollecnorm(PCM,DATA,[PAR,VAL])
% 	Normalize DATA according to PCM (Profile Classification Model)
%
% REQUIRED INPUTS:
% 	PCM: A Profile Classification Model structure
% 	DATA(Nz,Np): A collection of profiles
% 
% OPTIONAL INPUTS as [PAR,VAL]: List of optional parameters/values: 
% 	'usemodelnorm': Determine which data to use during the normalization step
% 		true,1 (default):  Use the model data
% 		false,0: Use DATA to apply the normalization method
% 
% Usage 2:
%
% DATAn = pcollecnorm(METHOD,DATA)
% 	Normalize DATA according to the METHOD
%
% REQUIRED INPUTS:
%	METHOD: An integer to choose the normalization method: 
% 		0: Nothing is done
% 		1 (default): Centered and standardized at each depth levels
% 		2: Centered at each depth levels
% 	DATA(Nz,Np): A collection of profiles
%
%
% OUTPUTS:
%	DATAn(Nz,Np): Normalized data
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

%- Default parameters:
usemodelnorm = false;

%- Load user parameters:
rnin = 2; % Number of required inputs
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

%- DATA size
[Nz,Np] = size(Xi);

%- Determine if we work with a PCM or not
if isa(ARG1,'struct')
	MODEL  = ARG1;
	METHOD = MODEL.normalization;
	if Nz ~= MODEL.DPTmodel
		error('Wrong DATA depth levels with regard to the PCM depth axis')
	end% if 
else
	METHOD = ARG1;
	usemodelnorm = 0;
end% if 

%- Normalize this !
switch METHOD
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

end %functionpcollecnorm