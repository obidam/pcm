function NewMatrix = covarrotate(Matrix, theta, ScaleRatio)
% covarrotate Rotate a squared covariance matrix
%
% rC = covarrotate(C, THETA) Rotate a covariance matrix
%
% Inputs:
%	C(D,D): Squared covariance matrix to rotate
% 	THETA: The angle in deg. to rotate (positive, clockwise)
% 
% Outputs:
%	rC: The C matrix rotated by THETA
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017,  OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-03-11 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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

M1 = [cos(theta) -sin(theta); sin(theta) cos(theta)];
NewMatrix = ScaleRatio^0 * M1 * Matrix * M1';

end %functioncovarrotate