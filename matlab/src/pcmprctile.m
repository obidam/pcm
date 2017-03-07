function varargout = pcmprctile(MODEL,DPT,DATA,M,varargin)
% pcmprctile Compute the percentile at each depth levels from a PCM and a given collection of profiles
%
% Y = pcmprctile(PCM,DPT,DATA,M,[PAR,VAL])
%
% REQUIRED INPUTS:
%	PCM (struct): A Profile Classification Model structure
% 	DPT(Nz,1): Depth axis of DATA (used to call pcmpredict and get the activation values)
% 	DATA(Nz,Np): A collection of profiles
% 	M: scalar or vector of percent values between 0 and 100
% 
% OPTIONAL INPUTS as [PAR,VAL]: List of optional parameters/values: 
% 	REDUCED: 0 (default) / 1: Indicate if DATA is reduced or not. If it's reduced
% 		real data will be recomposed using MODEL info.
% 	WEIGHT(Np,K): Possibly provide the profile weights to use to compute histograms.
% 		If not provided, DATA is classified using pcmpredict and ACTIVATIONS are used.
% 	CLASS: A list of integers to specify the class to determine
%
% OUTPUTS:
% 	Y: Percentiles of the values in DATA and CLASS
%
% EG:
% Ymedian = pcmprctile(PCM,DPT,DATA,50);
% plot(Ymedian,DPT);
% 
% Ymedian = pcmprctile(PCM,DPT,DATA,50,'CLASS',2);
% plot(Ymedian,DPT);
% 
% See Also: pcmpdfz
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-04-15 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)
% Revised: 2016-04-25 (G. Maze) Fixed a bug for the case where no data are defined at a given depth

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
REDUCED = false;

%- Load user parameters:
rnin = 4;
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

%- Get the correct number of vertical levels:
switch REDUCED
	case 0
		[Nz,Np] = size(DATA);
	case 1
		Nz = length(MODEL.DPTmodel);
end% switch 

%- Get the list of CLASS to diagnose:
if ~exist('CLASS','var')
	Klist = 1:MODEL.K;
else
	Klist = CLASS;
	if ~isempty(find(Klist>MODEL.K))
		error('CLASS numbers must be lower or equal to the PCM class number !');
	end% if 
	if ~isempty(find(Klist<0))
		error('CLASS numbers must be strictly positive !');
	end% if 
end% if 

%- Get profile weighting coefficients
if ~exist('WEIGHT','var')
	% Use activations as weights
	if ~REDUCED
		[~,~,WEIGHT] = pcmpredict(MODEL,DPT,DATA);
	else
		WEIGHT = gmmactiv(MODEL.mix,DATA');
	end% if
else
	% WEIGHT are provided as inputs
	if size(WEIGHT,1) ~= size(DATA,2) | size(WEIGHT,2) ~= MODEL.K
		error('WEIGHT must (Np,K)')
	end% if 
end% if 

%- Compute percentile at each depth levels
% Using activation-weighted statistics
Y = zeros(Nz,length(M),length(Klist));

for iz = 1 : Nz
	%-- Check if we need to handle reduced data or not
	switch REDUCED
		case 0 %--- Nope
			ds = DATA(iz,:);			
		case 1 %--- Yes
			% Recompose data for this depth level from EOFs		
			ds = MODEL.EOFs(iz,:)*DATA + MODEL.X_ref(iz);	
	end% switch

	%-- Compute stats for each class:
	for ik = 1 : length(Klist)
		iok = find(~isnan(ds(:)) & ~isnan(WEIGHT(:,Klist(ik)))); % Simpler trick to handle NaN in DATA or WEIGHT
		if ~isempty(iok)
			for ip = 1 : length(M)
				Y(iz,ip,ik) = wprctile(ds(iok),M(ip),WEIGHT(iok,Klist(ik)));
			end% for ip
		else
			Y(iz,1:length(M),ik) = NaN.*ones(1,length(M),1);
		end% if 
	end% for ik

end% for iz

%- Outputs
varargout(1) = {Y};

end %functionpcmprctile




















