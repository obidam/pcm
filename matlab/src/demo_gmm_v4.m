function varargout = demo_gmm_v4(varargin)
% demo_gmm_v4 Another Demonstrate of GMM
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
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


%- Load Netlab toolbox:
addpath('~/matlab/Netlab3.3');

%- Misc
range = @(X) [min(X(:)) max(X(:))];
hline = @() line(get(gca,'xlim'),[0 0]);
vline = @() line([0 0],get(gca,'ylim'));

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
prior = [0.5 0.5];

% Define dummy clusters:
% Clearly separated clusters:
c  = [2.0, 3.0; 1.0 1.0];  % Cluster centres in the 2D space:
sd = [1 0.1; 0.1 1]; % Cluster standard deviations
sd = [.5 0.5; 0.1 .5]; % Cluster standard deviations

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
	em_gmm_v2(data',K, 'covartype',covar_type,'debug',0); % Data should be Nd / Np

%- Plot results:
figure;  hold on
axl = ceil(max(abs(data(:)))); % Range of data (for axis limits)
[~,idm] = max([diff(range(data(:, 1))) diff(range(data(:, 2)))]);
axl = range(data(:, idm));

%-- Plot dataset overlaid with class elliptic Gaussian pdf contours
plot(data(:, 1), data(:, 2), 'k.','markersize',10,'color',[1 1 1]/2)
set(gca,'xlim',axl,'ylim',axl)
grid on, box on, axis square,hline;vline;
%axis([-1 1 -1 1]*axl),
%set(gca,'xtick',-axl:axl,'ytick',-axl:axl);

% Plot class Gaussian ellipses main axis:
col = 'rb'; 
for k = 1 : K
%	[l p] = plot_GMMellipse(MIX,1,2,k);
%	set([l p],'color',col(k))
	
	% Class covariance matrix:
	C = MIX.covars([1 2],[1 2],k);

	% Covariance eigen vectors 
	[v,d] = eig(C);

	% Ensure that eigenvector has unit length and plot ellipse axis
	col = 'kk';
	for idir = 1 : 2	
		v(:,idir) = v(:,idir)/norm(v(:,idir));
		start = MIX.centres(k,[1 2])-sqrt(d(idir,idir))*(v(:,idir)');
		endpt = MIX.centres(k,[1 2])+sqrt(d(idir,idir))*(v(:,idir)');
		linex = [start(1) endpt(1)];
		liney = [start(2) endpt(2)];
		l(idir) = line(linex, liney, 'Color', col(idir), 'LineWidth', 2);
	end% for idir
	
	% Class covariance matrix:
	C = MIX.covars([1 2],[1 2],k);

	% Possibly rotate Gaussian axis:
	% Rotation preserve eigen values and variance amplitudes
	if 1
		theta = 20*pi/180;
		M1 = [cos(theta) -sin(theta); sin(theta) cos(theta)];
		C  = M1 * C * M1';
	end% if 
	
	% Covariance eigen vectors 
	[v,d] = eig(C);

	% Ensure that eigenvector has unit length and plot ellipse axis
	col = 'rr';
	for idir = 1 : 2	
		v(:,idir) = v(:,idir)/norm(v(:,idir));
		start = MIX.centres(k,[1 2])-sqrt(d(idir,idir))*(v(:,idir)');
		endpt = MIX.centres(k,[1 2])+sqrt(d(idir,idir))*(v(:,idir)');
		linex = [start(1) endpt(1)];
		liney = [start(2) endpt(2)];
		l(idir) = line(linex, liney, 'Color', col(idir), 'LineWidth', 2);
	end% for idir
	
end% for
title(sprintf('Trained GMM with a %s covariance matrix\nand K=%i',MIX.covar_type,MIX.ncentres),'fontweight','bold','fontsize',14)







