function MODEL = classpermute(MODEL,NEWORDER)
% classpermute Permute class order in a Profile Classification Model (PCM)
%
% PCM = classpermute(PCM,NEWORDER);
%
% REQUIRED INPUTS:
%	PCM: A Profile Classification Model structure
% 	NEWORDER: The vector of integers from 1 to K of the new class order
% 
% OUTPUTS:
%	PCM: The updated Profile Classification Model structure
% 
% EG: 
% Create a PCM where classes are ordered by increasing SST:
%	[~, NEWORDER] = rename_labels(LABELS,POST,SST,'ascend'); 
%	PCM_sst = classpermute(PCM,NEWORDER);
%
% See Also: rename_labels
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-04-15 (G. Maze, Ifremer, Laboratoire d'Oceanographie Physique et Spatiale)

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

%- Check inputs
Klist = NEWORDER;
if ~isempty(find(Klist>MODEL.K))
	error('CLASS numbers must be lower or equal to the PCM class number');
end% if 
if ~isempty(find(Klist<0))
	error('CLASS numbers must be strictly positive');
end% if
if length(Klist) ~= MODEL.K
	error('It is not allowed to reduce the number of classes')
end% if 

%- Re-order
MODEL.mix.priors  = MODEL.mix.priors(Klist);
MODEL.mix.centres = MODEL.mix.centres(Klist,:);
switch MODEL.mix.covar_type
	case 'full'
		MODEL.mix.covars = MODEL.mix.covars(:,:,Klist);
	case 'diag'
		MODEL.mix.covars = MODEL.mix.covars(Klist,:);
	case 'spherical'
		MODEL.mix.covars = MODEL.mix.covars(1,Klist);		
end% switch 

end %functionclasspermute