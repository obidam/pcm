function [l p] = plot_GMMellipse(GMM,idref,id,k)
% plot_GMMellipse Plot Gaussian mode ellipse for a given GMM
%
% [ax el] = plot_GMMellipse(MIX,IDREF,ID,K)
%
% Inputs:
% 	MIX: The GMM model
% 	IDREF: The dimension to use along the x-axis
% 	ID: The dimension to use along the y-axis
% 	K: The Gaussian mode to plot
%
% Outputs:
% 	ax: Ellipse axes line handles
% 	el: Ellipse contours plot handles
%
% PCM is Profile Classification Modelling
% Copyright (C) 2014-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2014-12-10 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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

switch GMM.covar_type
	case 'full',      [v,d] = eig(GMM.covars([idref id],[idref id],k));
	case {'diag','mix'}, [v,d] = eig(diag(GMM.covars(k,[idref id])));
	case 'spherical', [v,d] = eig(diag([1 1]*GMM.covars(k)));
end% switch
%stophere

ALPH = 'ABCDEFG';
for idir = 1 : 2	
	% Ensure that eigenvector has unit length
	v(:,idir) = v(:,idir)/norm(v(:,idir));
	start = GMM.centres(k,[idref id])-sqrt(d(idir,idir))*(v(:,idir)');
	endpt = GMM.centres(k,[idref id])+sqrt(d(idir,idir))*(v(:,idir)');
	if 1
		col = 'kk';
		linex = [start(1) endpt(1)];
		liney = [start(2) endpt(2)];
		l(idir) = line(linex, liney, 'Color', col(idir), 'LineWidth', 2);
	else
		col = 'rb';
		linex = [start(1) endpt(1)];
		liney = [start(2) endpt(2)];
		l(idir) = line(linex, liney, 'Color', col(idir), 'LineWidth', 2);
		text(GMM.centres(k,idref),GMM.centres(k,id),'C','Color', 'k','fontsize',14);
		text(start(1),start(2),[ALPH(idir) '-1'],'Color', col(idir),'fontsize',14);
		text(endpt(1),endpt(2),[ALPH(idir) '-2'],'Color', col(idir),'fontsize',14);
	end% if 
end

% Plot ellipses of one standard deviation
theta = 0:0.02:2*pi;
x = sqrt(d(1,1))*cos(theta);
y = sqrt(d(2,2))*sin(theta);
% Rotate ellipse axes
ellipse = (v*([x; y]))';
ellipse2 = (v*([x; y]))';
% Adjust centre
ellipse  = ellipse + ones(length(theta), 1)*GMM.centres(k,[idref id]);
ellipse2 = 2*ellipse2 + ones(length(theta), 1)*GMM.centres(k,[idref id]);
p(1) = plot(ellipse(:,1), ellipse(:,2), 'k-','LineWidth', 2);
p(2) = plot(ellipse2(:,1), ellipse2(:,2), 'k-','LineWidth', 1);

%in1s = inpolygon(OBS(idref,labels==k),OBS(id,labels==k),ellipse(:,1), ellipse(:,2));
%in2s = inpolygon(OBS(idref,labels==k),OBS(id,labels==k),ellipse2(:,1), ellipse2(:,2));
%leg(1) = {sprintf('1 std (%0.1f%%)',length(find(in1s==1))*100/Klen(ik))};
%leg(2) = {sprintf('2 std (%0.1f%%)',length(find(in2s==1))*100/Klen(ik))};
%legend(p,leg,'location','best');
stophere

end %functionplot_GMMellipse