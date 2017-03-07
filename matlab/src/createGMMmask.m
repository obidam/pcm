function MASK = createGMMmask(LSmask,FIELD,DPT,varargin)
% createGMMmask Create a mask from land/sea and 3D field data constraints
%
% MASK = createGMMmask(LSmask,FIELD,DPT) 
% This function will create a mask with:
% - 1 where LSmask is 1 & FIELD is 1 at first level
% - 0 elsewehere
%
% MASK = createGMMmask(LSmask,FIELD,DPT,'Zmin',VALUE) 
% This function will create a mask with:
% - 1 where LSmask is 1 & FIELD is 1 at all levels above minimum depth VALUE
% - 0 elsewehere
%
% Inputs:
% 	LSmask (lat/lon): A lat/lon mask field with 0/1 values
% 	FIELD: The dpt/lat/lon 3D field to use for 0/1 extraction
% 	DPT: The vertical axis (negative, downward) of FIELD
% Options:
% 	Zmin: A minimum depth for FIELD values to be defined (0 by default)
%
% Outputs:
% 	MASK: A (lat/lon) mask field of 0/1 values with 1 where data should be kept for classification 
%
% PCM is Profile Classification Modelling
% Copyright (C) 2015-2017,  OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2015-11-18 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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

%- Init
Zmin = 0;

%- Load user options:
if nargin-3 > 0
	if mod(nargin-3,2) ~= 0
		error('Options must come in pairs: OPTION,VALUE')
	end% if 
	for in = 1 : 2 : nargin-3
		eval(sprintf('%s = varargin{in+1};',varargin{in}));
	end% for in	
	clear in
end% if

%- More init
[nz ny nx] = size(FIELD);
d  = permute(FIELD,[1 3 2]);
da = d(:,1:nx*ny);

%- Apply LSmask:
c = LSmask';
ikeep1 = c(:)==1;

%- Then minimum depth:
Zi = repmat(DPT(:),[1 nx*ny]);
Zi(isnan(d(:,1:nx*ny))) = NaN;
Zi = min(Zi)';
ikeep2 = Zi <= Zmin; 

%- Create lat/lon mask
ii = zeros(1,nx*ny);
ii(ikeep1 & ikeep2) = 1;
MASK = zeros(nx,ny)*NaN;
MASK(ikeep1&ikeep2) = 1;
MASK = MASK'; % Ny,Nx

end %functioncreateGMMmask