function varargout = pcmpdfz(MODEL,DPT,DATA,M,varargin)
% pcmpdfz Compute the PDF at each depth levels (PDFz) from a PCM and a given collection of profiles
%
% [PDFz,X,PDFzO,PDFzM] = pcmpdfz(PCM,DPT,DATA,M,[PAR,VAL])
%
% REQUIRED INPUTS:
%	MODEL (struct): A Profile Classification Model structure
% 	DPT(Nz,1): Depth axis of DATA (used to call pcmpredict and get the activation values)
% 	DATA(Nz,Np): A collection of profiles
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
% 	[PDFz,X,PDFzO,PDFzM] = pcmpdfz(MODEL,DPT,DATA,range(DATA(:),100));
% 	pcolor(X,DPT,squeeze(PDFz(:,:,3))); title('Class 3 PDFz')
% 
% See Also: pcmprctile
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

%- Define PDF bin axis
if length(M) == 1
	if ~REDUCED
		rg = range(DATA(:),M);
	else
		error('I don''t know yet how to do determined a range with reduced data !')
	end% if 
else
	rg = M;
end% if 

dr = unique(diff(rg));
if length(dr) > 1 & ~isempty(find(diff(dr)>1e-14))
	% Check if difference are significant and not epsilonesc
	error('I don''t know yet how to do this with an irregular set of bins !')
else
	dr = dr(1);
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

%- Compute PDF at each depth levels
% Using activation-weighted statistics
PDFz  = zeros(Nz,length(rg),length(Klist));
PDFzO = zeros(Nz,length(rg));
PDFzM = zeros(Nz,length(rg));

for iz = 1 : Nz
	%-- Check if we need to handle reduced data or not
	switch REDUCED
		case 0 %--- Nope
			ds = DATA(iz,:);			
		case 1 %--- Yes
			% Recompose data for this depth level from EOFs		
			ds = MODEL.EOFs(iz,:)*DATA + MODEL.X_ref(iz);	
	end% switch

	%-- Compute stats for the entire collection:
	c = histc(ds,rg); 
	c = c/sum(c*dr); % Transform weighted histogram into a pdf (integral is one, not nb of points)
	PDFzO(iz,:) = c;

	%-- Compute stats for each class:
	for ik = 1 : length(Klist)
		c = histcw(ds,WEIGHT(:,Klist(ik)),rg); % Count data weighted by WEIGHTS
		c = c/sum(c*dr); % Transform weighted histogram into a pdf (integral is one, not nb of points)
		PDFz(iz,:,ik) = c;
		PDFzM(iz,:) = PDFzM(iz,:) + MODEL.mix.priors(Klist(ik))*squeeze(PDFz(iz,:,ik));
	end% for ik

end% for iz

%- Outputs
varargout(1) = {PDFz};
varargout(2) = {rg};
varargout(3) = {PDFzO};
varargout(4) = {PDFzM};

end %functionpcmpdfz


% Compute a weighted histogram:
function [histw histv] = histcw(v, w, edges) 
	nx  = length(edges);
	iin = find(v>=edges(1) & v<=edges(end)); % To avoid any subs == 0 that would create errors in accumarray
	X = v(iin);
	W = w(iin);
	[N, subs] = histc(X,edges); % N(i) = sum(subs==i); subs is zero for out of range values
	histv = accumarray(subs(:),ones(size(subs(:))),[nx,1]); % Classic histogram
	histw = accumarray(subs(:),W,[nx,1]); % Weighted histogram
end%end function


































