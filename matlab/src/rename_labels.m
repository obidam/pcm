function [newlabels knew wm0] = rename_labels(labels,weights,scalar,mode)
% reorder_labels Re-name labels IDs
%
% [NEW_LABELS KLIST VAL] = rename_labels(LABELS,WEIGHT,SCALAR,MODE) 
%
% Create new labels according to the WEIGHT-weighted mean value of SCALAR
% of each class. MODE selects the direction of the sort.
% 
% Inputs:
%	LABELS(Np,1): Labels of the profiles, an integer between 1 and K
% 	WEIGHT(Np,1): A weight vector
% 	SCALAR(Np,1): A scalar value for each profile
% 	MODE: A string indicating the direction of the sort:
% 		'ascend' results in ascending order
%		'descend' results in descending order
% 
% Outputs:
%	NEW_LABELS(Np,1): The new labels
% 	KLIST(1,K): The old to new Class id ordering
% 	VAL(1,K): The K WEIGHT-weighted mean values of SCALAR used to re-order labels.
% 
% Eg:
% 	Create new labels so that the lowest class label has the smallest mean SST:
% 		[oLABELS Klist Ksst] = rename_labels(LABELS,ACTI,SST,'ascend');
% 	
% See Also: classpermute 
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2015-01-25 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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


% Determine the number of class:
K = length(unique(labels(~isnan(labels))));

% Compute for each class the weighted mean value
for ik = 1 : K
	ip = find(labels==ik); % Profiles attributed to this class
	wm0(ik) = wmean(weights(ip),scalar(ip));
end% for ik

% Sort weighted mean values
[wm knew] = sort(wm0,mode);

% Create new label vector
newlabels = zeros(size(labels))*NaN;
for ik = 1 : K
	ip = find(labels==knew(ik));
	newlabels(ip) = ik;
end% for ik

end %functionreorder_labels