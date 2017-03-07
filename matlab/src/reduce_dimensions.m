function varargout = reduce_dimensions(X,THRESHOLD,varargin)
% reduce_dimensions Reduce the dimension of a dataset with PCA/SVD
%
% [Xr,B,V,L,Xref,Xp,Xrnorm,Breal] = reduce_dimensions(X,EXPVAR)
% 
% Reduce the dimension of a dataset with PCA. 
% X is a collection of profiles with a large number of vertical 
% levels. This routine is used to reduced that number.
%
% Inputs:
% 	X [Nz,Np] : a collection of Np profiles (cols) with Nz vertical levels (rows).
% 	EXPVAR : the fraction of variance between 0 and 100% that must be 
% 		retained in the reduced dataset.
%
% Outputs:
%	Xr [Nc,Np] : the collection of Np profiles in the reduced Nc-dimensional 
% 		space of EOFs, ie the Principal component scores or PCs.
% 	B [Nz,Nc] : profiles defining the new Nc-dimensional space, ie the 
% 		Principal component coefficients or EOFs.
% 	V [1,Nc] : the fraction of variance explained by individual axis of the
% 		new space.
% 	L [1,Nc] : Principal component variances, ie eigenvalues of the covariance 
% 		matrix.
% 	Xref [Nz,1] : is the sample mean of the input collection X, used as as reference
% 		to center the collection.
%  	Xp [Nz,Np] : the re-composed collection, to be compared with X.
%	Xrnorm [Nc,Np] : Normalized Principal component scores.
% 	Breal [Nz,Nc] : profiles defining the new Nc-dimensional space with their 
% 		dimensional units.
%
% The re-composed collection is the backward projection of the reduced dataset onto
% the original Nz-dimensional space: Xp = B*Xr + repmat(Xref,[1 size(Xr,2)]);
%
% PCM is Profile Classification Modelling
% Copyright (C) 2015-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2015-01-21 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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

% We use notation from the Hartmann, 2002.
% EOF Analysis via Singular Vector Decomposition of the Data Matrix, page 61.
% See also the section "EOF Analysis Based on Normalized Data" 

% Assume that X is of size: structure/space vs sampling/profile
[Nz, Np] = size(X); % So that X is (Nz,Np)

% Max number of EOFs to compute:
if isinf(THRESHOLD)
	Nmax = Nz;
else
	Nmax = 20;
end% if 

% Load user options:
if nargin-2 > 0
	if mod(nargin-2,2) ~= 0
		error('Options must come in pairs: OPTION,VALUE')
	end% if 
	for in = 1 : 2 : nargin-2
		eval(sprintf('%s = varargin{in+1};',varargin{in}));
	end% for in	
	clear in
end% if
if Nmax > Nz
	warning(sprintf('The number of eofs (%i) exceeds the possible limit (%i)',Nz,Nmax))
	Nmax = Nz;
end% if 

%indexes_training = 1:Np;
%indexes_training = unique(randi(size(X,2),[1,floor(size(X,2)/2)]));
%indexes_test = setdiff(indexes,indexes_training);
%X = matrix(indexes_training,:)'; Np = length(indexes_training);
X_samplemean = mean(X,2);
if size(X,1) > size(X,2) % This is not what we expect: Nz > Np
	warning('You should verify the size of data matrix, it should be (Nz,Np)')
end% if 
if ~isempty(find(isnan(X_samplemean)==1))
	error('Data cannot contain NaN to be reduced')
end% if 

% Remove the sample mean profile:
Xtilde = X - repmat(X_samplemean,[1 Np]);
	
% Find eigenvectors and singular values of normalized data Xtilde, such as:
% Xtilde = U * S * V';
% - U contains eigenvectors of the structural covariance matrix X*X', ie taking the inner 
% product over sampling/profile leaving the covariance between structure/space points.
% - V contains eigenvectors of the sample covariance matrix X'*X, ie taking the inner 
% product over structure/space leaving the covariance between sampling/profile points.
% - S (singular values) contains sqrt(eigenvalues) of covariance X*X' and X'*X.
% try
	[U,S,V] = svds(Xtilde,Nmax); % Retains the n largest components
% catch ME
% 	stophere
% end% try

% Compute the Principal Components (Nc,Np) by projecting structural EOFs on data:
Z = U'*Xtilde;  % (Nc,Np)
Ztilde = S^(-1/2)*Z; % Normalize PCs

% Compute EOFs with their dimensional units by projection of the normalized
% PCs on the original data.
D = X*Ztilde'/Np; % (Nz,Nc)

% Compute eigenvalues from singular values:
% Eigenvalues contain the variance explained by each component
lambda = S^2;

% Compute fraction of explained variance by each component
% The fraction of the total variance explained by a particular eigenvector 
% is equal to the ratio of that eigenvalue to the trace of the eigenvalue 
% matrix, which is equal to the trace of the covariance matrix.
explained = diag(lambda)/trace(lambda); % A fraction between 0 and 1

% Link to Matlab pca routine vocabulary:
if 1
	coeff  = U;   % Principal component coefficients (EOFs, eigenvectors)
	score  = Z';  % Principal component scores (PCs)
	coeff_scaled = D;
	score_scaled = Ztilde';
else
	coeff  = D;        % Principal component coefficients (EOFs, eigenvectors) with their dimensional units
	score  = Ztilde';  % Normalized Principal component scores (PCs)
end% if 
latent = diag(lambda)'; % Principal component variances (that is the eigenvalues of the covariance matrix)

%Ztest = U'*(matrix(indexes_test,:)'-repmat(X_samplemean,[1 length(indexes_test)]));

% Identify the max number of EOFs to keep:
if isinf(THRESHOLD)
	important_components = Nmax;
else
	important_components = find(((cumsum(explained)))>=THRESHOLD/100, 1, 'first');	
	if important_components == Nmax
		error(sprintf('The variance level of %0.4f cannot be reached with the default number of EOFs (%i)\nTry to reduce the variance level or use Inf to compute the complete set of PCA.',THRESHOLD,Nmax));
	end% if 
end% if

% and squeeze EOFs and PCs to important_components:
projection = score(:,1:important_components)';
profiles   = (projection'*coeff(:,1:important_components)'+repmat(X_samplemean,[1 Np])')';

%- Output
%	Xr [Nc,Np] : the collection of Np profiles in the reduced Nc-dimensional 
% 		space of EOFs, ie the Principal component scores or PCs.
varargout(1) = {projection};

% 	B [Nz,Nc] : profiles defining the new Nc-dimensional space, ie the 
% 		Principal component coefficients or EOFs.
varargout(2) = {coeff(:,1:important_components)};

% 	V [1,Nc] : the fraction of variance explained by individual axis of the
% 		new space.
varargout(3) = {explained(1:important_components)'};

% 	L [1,Nc] : Principal component variances, ie eigenvalues of the covariance 
% 		matrix.
varargout(4) = {latent(1:important_components)};

% 	Xref [Nz,1] : is the sample mean of the input collection X, used as as reference
% 		to center the collection.
varargout(5) = {X_samplemean};

%  	Xf [Nz,Np] : the re-composed collection, to be compared with X.
varargout(6) = {profiles};

%	Xrnorm [Nc,Np] : the collection of Np profiles in the reduced Nc-dimensional 
% 		space of EOFs, ie the Normalized Principal component scores or PCs.
varargout(7) = {score_scaled(:,1:important_components)'};

% 	Breal [Nz,Nc] : profiles defining the new Nc-dimensional space, ie the 
% 		Principal component coefficients or EOFs with their dimensional units
varargout(8) = {coeff_scaled(:,1:important_components)};

end %functionreduce_dimensions