function varargout = demo_gmm_v3(varargin)
% demo_gmm_v3 Same as demo_gmm_v2 but with interactive nb of cluster
% Note that you need the TEOS-10 library to run this demo !
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

%- Load Netlab toolbox:
addpath('~/matlab/Netlab3.3');
try
	addpath(genpath('~/matlab/gsw')) % New TEOS-10 library (gsw_*) to replace seawater
catch
	error('You need the TEOS-10 library to run this demo !');
end

%- Parameters
Kmin = 1; 
Kmax = 10;
K = 4;
MIX = cell(1,Kmax-Kmin+1);

%- Load data:
% We create a dataset with T/S at z=-300 in the North Atlantic with Argo profils
load('demo_gmm_v2_data.mat');
iz = find(S.dpt<=-300,1); % 
data(1,:) = S.psal(iz,:);
data(2,:) = S.temp(iz,:);

[Nd Np] = size(data); % Nb od dimensiosn / Nb of profils

%- Pre-process data:
% Convert to TEOS10: 
pres = gsw_p_from_z(-abs(S.dpt(iz))*ones(size(S.lat)),S.lat); % Pressure
data(1,:) = gsw_SA_from_SP(data(1,:),pres,S.lon,S.lat); % Absolute Salinity
data(2,:) = gsw_CT_from_t(data(1,:),data(2,:),pres); % Conservative Temperature

% Scale to density:
[~, alpha, beta] = gsw_rho_alpha_beta(data(1,:),data(2,:),pres);
data(1,:) = data(1,:).*beta;
data(2,:) = data(2,:).*alpha;

% Standardize/center
Xscl_ave = mean(data,2); 
Xscl_std = std(data,[],2);
data = (data-repmat(Xscl_ave,[1 Np]))./repmat(Xscl_std,[1,Np]);

%- Train a GMM on the dataset:
MIX{K} = trainthis(data,K);

%- Plot it
figure; hold on
% For faster plot rendering, we only use a fraction of the dataset
ip = randi(Np,1,1000);
axl = ceil(max(abs(data(:)))); % Range of data (for axis limits)

cf = get(0,'currentfigure');
setappdata(cf,'gmm',MIX);
setappdata(cf,'gmmdata',data);

plot(data(1,ip), data(2,ip), 'k.','markersize',10,'color',[1 1 1]/2)
grid on, box on, axis square,axis([-1 1 -1 1]*axl),
set(gca,'xtick',-axl:axl,'ytick',-axl:axl);
plotthisgmm(MIX{K});
t = text(-axl(1),axl(1),sprintf('llh=%0.5f',MIX{K}.llh),'tag','LLH');
set(t,'VerticalAlignment','top','HorizontalAlignment','left','fontsize',14)

title(sprintf('Use slider to change the nb of components in GMM'),'fontweight','bold','fontsize',14)

%- Add slider
jSlider = javax.swing.JSlider; [jhSlider, hContainer] = javacomponent(jSlider,'East');
set(jSlider,'Orientation',jSlider.VERTICAL); 
set(jSlider,'Minimum',Kmin,'Maximum',Kmax,'Value',K);

set(jSlider,'paintTicks',true)
set(jSlider,'PaintLabels',true);
set(jSlider,'snapToTicks',1)
set(jSlider,'majorTickSpacing',1);
set(jSlider,'minorTickSpacing',1);

%- Set Callback:
cbk = @(obj,eve) cbk2(obj,eve);
%cbk = @(obj,eve) set(findobj(ax,'tag','interactive'),'ydata',fct(y,get(obj,'value')));
hjSlider = handle(jSlider, 'CallbackProperties');
set(hjSlider, 'StateChangedCallback', cbk);  %alternative

end %functiondemo_gmm_v3


function cbk2(obj,eve)

	cf = get(0,'currentfigure'); % this should be determined from obi	
	MIX  = getappdata(cf,'gmm');
	data = getappdata(cf,'gmmdata');
	k = fix(get(obj,'value'));
	if isempty(MIX{k})
		%disp('Updating GMM, please wait before chaging slider again ...')
		MIX{k} = trainthis(data,k);
		setappdata(cf,'gmm',MIX);		
	end% if 
	
	plotthisgmm(MIX{k});
	set(findall(cf,'tag','LLH'),'string',sprintf('llh=%0.3f',MIX{k}.llh));
	%disp('Done')
	
end

function MIX = trainthis(data,K)
	[POST PROB ACTI LLH LABELS MIX ERRLOG MIXKMEAN ERRLOGKMEAN] = ...
		em_gmm_v2(data,K,'covartype','full','debug',0); % Data should be Nd / Np	
	MIX.llh = LLH;
end

function h = plotthisgmm(MIX)

	delete(findall(get(0,'currentfigure'),'tag','interactive'));	
	K = MIX.ncentres;
	col = hsv(K);
	for k = 1 : K
		[l p] = plot_GMMellipse(MIX,1,2,k);
		set([l p],'color',col(k,:),'tag','interactive')
	end% for
	h = [l p];
	
end% function

