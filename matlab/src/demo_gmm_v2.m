function varargout = demo_gmm_v2(varargin)
% demo_gmm_v2 Demonstrate GMM in 2D with a T/S diagram
% Note that you need the TEOS-10 library to run this demo !
% 
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-03-08 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)
% Revised: 2017-01-18 (G. Maze) Fixed a bug with em_gmm_v2 call

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

clear

%- Load toolbox:
addpath('~/matlab/Netlab3.3');
try
	addpath(genpath('~/matlab/gsw')) % New TEOS-10 library (gsw_*) to replace seawater
catch
	error('You need the TEOS-10 library to run this demo !');
end	

%- Load data:
% We create a dataset with T/S at z=-300 in the North Atlantic with Argo profils
load('demo_gmm_v2_data.mat');
iz = find(S.dpt<=-300,1); % 
data(1,:) = S.psal(iz,:);
data(2,:) = S.temp(iz,:);

[Nd Np] = size(data); % Nb od dimensiosn / Nb of profils

% Convert to TEOS10: 
pres = gsw_p_from_z(-abs(S.dpt(iz))*ones(size(S.lat)),S.lat); % Pressure
data(1,:) = gsw_SA_from_SP(data(1,:),pres,S.lon,S.lat); % Absolute Salinity
data(2,:) = gsw_CT_from_t(data(1,:),data(2,:),pres); % Conservative Temperature

data0 = data; % Store original data for final plots

%- Pre-process data:
% Scale to density:
if 0
	[~, alpha, beta] = gsw_rho_alpha_beta(data(1,:),data(2,:),pres);
	data(1,:) = data(1,:).*beta;
	data(2,:) = data(2,:).*alpha;
end% if 

% Standardize/center
Xscl_ave = mean(data,2); 
Xscl_std = std(data,[],2);
data = (data-repmat(Xscl_ave,[1 Np]))./repmat(Xscl_std,[1,Np]);

%- Train a GMM on the dataset:
covar_type = 'full'; % Co-variance matrix form, any in: 'spherical','diag','full'
K = 5;
[POST PROB ACTI LLH LABELS MIX ERRLOG MIXKMEAN ERRLOGKMEAN] = ...
	em_gmm_v2(data,K,'covartype',covar_type,'debug',0); % Data should be Nd / Np


%- Plot results:
% For faster plot rendering, we only use a fraction of the dataset
ip = randi(Np,1,1000);

figure; jw=2;iw=1+ceil(K/2);ipl=0;
axl = max(abs(data(:))); % Range of data (for axis limits)

%-- Plot dataset overlaid with class elliptic Gaussian pdf contours
ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
plot(data(1,ip), data(2,ip), 'k.','markersize',10,'color',[1 1 1]/2)
grid on, box on, axis square,axis([-1 1 -1 1]*axl),
%set(gca,'xtick',linspace(-axl,axl,20),'ytick',-axl:axl);

col = hsv(K);
for k = 1 : K
	[l p] = plot_GMMellipse(MIX,1,2,k);
	set([l p],'color',col(k,:))
end% for
title(sprintf('Trained GMM with a %s covariance matrix\nand K=%i',MIX.covar_type,MIX.ncentres),'fontweight','bold','fontsize',14)

%-- Plot Hard classification
% i.e. the class a point is attributed to.
ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
scatter(data(1,ip), data(2,ip), 20, LABELS(ip));
grid on, box on, axis square,axis([-1 1 -1 1]*axl),
%set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
caxis([0 K]); colorbar
title(sprintf('Data labels\nhard-classification'),'fontweight','bold','fontsize',14)

%-- Plot posteriors for each class:
% i.e. the probability for a datapoint to belong the class k
for k = 1 : K
	ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
	scatter(data(1,ip), data(2,ip), 20, POST(ip, k));
	grid on, box on, axis square,axis([-1 1 -1 1]*axl),
%	set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
	caxis([0 1]); colorbar
	title(sprintf('Posteriors for class %i: p(k=%i|c)\nPrior = %0.2f',k,k,MIX.priors(k)),'fontweight','bold','fontsize',14)
end% for k

%- Now convert GMM results in the real dimensional plan:
figure
ip = randi(Np,1,2000);
scatter(data0(1,ip), data0(2,ip), 20, LABELS(ip));
grid on, box on, axis square
caxis([0 K]); colorbar
title(sprintf('Data labels\nhard-classification'),'fontweight','bold','fontsize',14)

end %functiondemo_gmm_v2

	