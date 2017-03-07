function varargout = demo_gmm_v1(varargin)
% demo_gmm_v1 Simplest demonstrate of GMM in 2D
%
% Simplest demonstration ever of a GMM on a dummy dataset
% 
% PCM is Profile Classification Modelling
% Copyright (C) 2016,  OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-03-08 (G. Maze, Ifremer, Laboratoire d'Oceanographie Physique et Spatiale)

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

%- Load Netlab toolbox:
addpath('~/matlab/Netlab3.3');

%- Generate a dummy 2D dataset:
% mixture of 2 Gaussians in two dimensional space
Nd = 2; % Nb of dimensions

% Fix seed for reproducible results
randn('state', 42);

% Number of points in the dataset:
Np = 500;
data = randn(Np, Nd);

% Generation it:
K = 2; % Nb of cluster/class

% Priors for the ncentres clusters
% Prior is the probability of the cluster for the dataset, 
% i.e. the expected fraction of points attributed to the 
% class vs the dataset total nb of points
prior = [0.3 0.7];

% Define dummy clusters:
% Clearly separated clusters:
c  = [2.0, 3.0; 1.0 0.0];  % Cluster centres in the 2D space:
sd = [.3 .7; 0.5 0.5]; % Cluster standard deviations

% Overlapping clusters:
c  = [2.0, 2.0; 0.0 0.0];  % Cluster centres in the 2D space:
sd = [.4 .7; 0.4 0.1]; % Cluster standard deviations

% Identify points for each center
ind{1} = 1:prior(1)*Np;
for ic = 2 : K
	ip = 1:prior(ic)*Np;
	ind{ic} = ip+max(ind{ic-1});
end% for ic

for ic = 1 : K
	data(ind{ic},1) = data(ind{ic},1)*sd(ic,1) + c(ic,1);
	data(ind{ic},2) = data(ind{ic},2)*sd(ic,2) + c(ic,2);
end% for ic

% Clear workspace
clear ic ind

%- Train a GMM on the dataset:
covar_type = 'full'; % Co-variance matrix form, any in: 'spherical','diag','full'
[POST PROB ACTI LLH LABELS MIX ERRLOG MIXKMEAN ERRLOGKMEAN] = ...
	em_gmm_v2(data',K,'covartype',covar_type,'debug',0); % Data should be Nd / Np

%- Plot results:
figure; iw=2;jw=2;ipl=0;
axl = ceil(max(abs(data(:)))); % Range of data (for axis limits)

%-- Plot dataset overlaid with class elliptic Gaussian pdf contours
ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
plot(data(:, 1), data(:, 2), 'k.','markersize',10,'color',[1 1 1]/2)
grid on, box on, axis square,axis([-1 1 -1 1]*axl),
set(gca,'xtick',-axl:axl,'ytick',-axl:axl);

col = 'rb'; 
for k = 1 : K
	[l p] = plot_GMMellipse(MIX,1,2,k);
	set([l p],'color',col(k))
end% for
title(sprintf('Trained GMM with a %s covariance matrix\nand K=%i',MIX.covar_type,MIX.ncentres),'fontweight','bold','fontsize',14)

%-- Plot posteriors for each class:
% i.e. the probability for a datapoint to belong the class k
for k = 1 : K
	ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
	scatter(data(:, 1), data(:, 2), 20, POST(:, k));
	grid on, box on, axis square,axis([-1 1 -1 1]*axl),
	set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
	caxis([0 1]); colorbar
	title(sprintf('Posteriors for class %i: p(k=%i|c)\nPrior = %0.2f',k,k,MIX.priors(k)),'fontweight','bold','fontsize',14)
end% for k

%-- Plot Hard classification
% i.e. the class a point is attributed to.

ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
scatter(data(:, 1), data(:, 2), 20, LABELS);
grid on, box on, axis square,axis([-1 1 -1 1]*axl),
set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
caxis([0 K]); colorbar
title(sprintf('Data labels\nhard-classification'),'fontweight','bold','fontsize',14)
















