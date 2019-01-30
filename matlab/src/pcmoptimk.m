function varargout = pcmoptimk(NIND,DATA,DPT,MODE,varargin)
% pcmoptimk Determine the optimum number of classes K (using BIC)
%
% Determine the optimum number of class K as the one minimizing an
% ensemble mean of BIC estimates for a range of K.
% Each BIC value is estimated using NIND independent profiles randomly
% peaked out of DATA.
% 
% [Kbest, BIClist] = pcmoptimk(NIND,DATA,DPT,'informed',PCM,[PAR,VAL])
% 	Use the known PCM to pre-process DATA.
% 
% [Kbest, BIClist] = pcmoptimk(NIND,DATA,DPT,'blind',[PAR,VAL]) 
% 	The entire data set DATA is new
%
% REQUIRED INPUTS:
% 	NIND (integer): The number of independent profiles to consider
% 	DATA(Nz,Np): The field to classify with a GMM, Np profiles with Nz depth levels.
% 	
% OPTIONAL INPUTS as [PAR,VAL]:
% 	NENS (Default: 50): The size of the ensemble (nb of training), must be smaller than Np
% 	TRNOPT (no default, cell of strings): Optionals inputs to be used when calling PCMTRAIN
% 	KLIST (Default: [1:30]) The list of number of class to explore.
% 	PCM (no default): The PCM to be used when MODE='rookie'
%
%
% Help:
% for a 150km decorrelation scale over the North Atlantic (Ninove et al, 2016: fig 9)
% NindepX = lldist([1 1]*40,[-70 -15])/1000/150 ~= 31
% NindepY = lldist([15 55],[1 1]*-45)/1000/150 ~= 30
% Np = NindepX*NindepY
%
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-05-20  (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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


%- Default parameters:
NrunsKM = 5;
KLIST = 1:30;
COVARTYPE = 'full';
NENS = 50;
TRNOPT = {};

%- Load user parameters:
switch lower(MODE)
	case 'informed'
		rnin = 5; % Number of required inputs
		instart = 2;
		if ~isa(varargin{1},'struct')
			error('Using MODE ''informed'', you must provide a PCM')
		else
			PCM = varargin{1};
		end% if
	case 'blind'
		rnin = 4; % Number of required inputs
		instart = 1;
end% switch 
if nargin > rnin
    if mod(nargin-rnin,2) ~= 0
        error('Parameters must come in pairs: PAR,VAL');
    end% if
    for in = instart : 2 : nargin-rnin
        eval(sprintf('%s = varargin{in+1};',varargin{in}));
    end% for in
    clear in
elseif nargin == rnin
	% OK
else
	error('Bad number of parameters');
end% if
clear rnin instart

if isempty(intersect(lower(MODE),{'informed','blind'}))
	error('Invalid mode');
end% if 

%- 
[Nz Np] = size(DATA);
ipreg = randperm(NIND*fix(size(DATA,2)/NIND));

%- Pre-process data according to the running MODE:
switch lower(MODE)
	case 'blind'
		% We train a PCM with K=1 in order to get the reduced data set
		% as fast as possible
		[~,~,~,DATAr] = pcmtrain(DATA,1,COVARTYPE,DPT,TRNOPT{:});
		
	case 'informed'
		% We use the known PCM to reduce the data set:
		DATAr = pcmreduce(PCM,DPT,DATA);
end% switch 

%- Compute an ensemble of NENS of BIC(KLIST)
ii = 0;
for irun = 1 : NENS
	for itest = 1 : length(KLIST)
		ii = ii + 1;
		[irun itest]
%		jvmwaitbar(NENS*length(KLIST),ii,'Computing BIC statistics...')
		k  = KLIST(itest);
		Xr = DATAr(:,ipreg(randperm(length(ipreg),NIND)));
		try			
			[post, prob, acti, nllh, labels, mix, errlog, ~, errlogKM] = ...
				em_gmm_v2(Xr,k,TRNOPT{:});
			[B(itest,irun) B1(itest,irun) B2(itest,irun) N(itest,irun)] = gmmbic(mix,nllh,Xr);
		catch
			B(itest,irun)  = NaN;
			B1(itest,irun) = NaN;
			B2(itest,irun) = NaN;
			N(itest,irun)  = NaN;
		end
	end% for itest
end% for irun

%- Determine the optimum K:
[~,ik] = min(nanmean(B,2));
Kbest = KLIST(ik);

if nargout == 0
	
	%- Plot BIC
	figure;
	p = ploterr(KLIST,nanmean(B,2),0,nanstd(B,[],2));
	set(p,'color','k','linewidth',2)
	grid on, box on,xlim(range(KLIST))
	vline(Kbest,'color','r','linewidth',2);
	set(gca,'fontsize',14);
	xlabel('Number of class','fontweight','bold','fontsize',14);
	ylabel('BIC','fontweight','bold','fontsize',14);
	title(sprintf('BIC Statistics for %i runs with %i independent profiles',NENS,NIND),'fontweight','bold','fontsize',16)

else
	varargout(1) = {Kbest};
	varargout(2) = {B};
	varargout(3) = {KLIST};
end% if 

end %functionpcmoptimk