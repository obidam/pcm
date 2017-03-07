%DEF Help for the Profile Classification Model (PCM) structure
%
% Version 1.0
% A Profile Classification Model (PCM) structure can have the following fields:
% 
% GLOBAL ATTRIBUTES
% 	readme:  
% 	Np: Number of profiles used when the model was trained
% 
% INTERPOLATION 	
% 	DPTmodel: The vertical depth axis of the model (negative, downward)
%
% NORMALIZATION
% 	normalization
% 	X_ave
% 	X_std
% 
% REDUCTION
% 	doREDUCE: Do we compress the data with PCA or not ?
% 	EOFs: Eigenvectors of the reduced space
% 	X_ref: Center of the array used in PCA
% 	V: Variance of each of the PCA dimensions
% 	maxvar: The effective variance retainted within the reduced space
% 
% CLASSIFICATION MODEL
% 	covarTYPE:
% 	K: The number of class in the GMM model
%	LLH: Negative Log likelihood of the GMM
%	mix: The Netlab GMM structure
% 
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-04-18 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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


help pcm


