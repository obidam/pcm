function [post prob acti log_likelihood main mixGM errlogGM mix0 errlogKM] = em_gmm_v2(matrix,number_classes,varargin)
% em_gmm_v2 (active) Train a GMM model with EM algorithm and KMean initialization.
%
% [POST PROB ACTI LLH LABELS MIX ERRLOG MIXKMEAN ERRLOGKMEAN] = em_gmm_v2(DATA,K,[PAR,VAL]) 
%
% Inputs:
%	DATA double array [Nd,Np]: the dataset, a collection of Np vectors of Nd dimensions.
% 	K: an integer defining the number of components/classes of the GMM
% 	PAR/VAL: Possible pairs of parameter/variable to customize the function behavior:
% 		- covartype: a string to determine the type of Netlab GMM covariance type to use. It can be:
% 			- 'spherical', for isotrope variance matrix (identity*sigma)
% 			- 'diag', for diagonal covariance (dimensions not correlated, 1 sigma per dimension)
% 			- 'full' (default), for full covariance matrix (correlated dimensions)
% 		- initializationMode: 'kmean' or 'random'
% 		- NrunsKM (default: 20): an integer to fix the number of KMeans runs among which to choose the GMM initial conditions
% 		- NiterKM (default: 200): The maximum nb of iterations to perform during each KMean initialization runs
% 		- NiterGM (default: 1000): The maximum nb of iterations to perform during the GMM run
% 		- check_covar (default: 0)
% 		- debug (default: 0): Print out convergence figures
% 
% Outputs:
% 	POST: Posterior probabilities
% 	PROB: Mixture probabilities
% 	ACTI: Component probabilities
% 	LLH: Log likelihood of the model
% 	LABELS: Class labels attributed to data
% 	MIX: The Netlab GMM structure
% 	ERRLOG: Log of the model LLH during training iterations
% 	MIXKMEAN: The Netlab K-mean structure of the initial model
% 	ERRLOGKMEAN: Log of the initial model LLH during training K-mean iterations
% 	
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2015-11-26 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)
% Revised: 2016-04-13 (G. Maze) Moved COVARTYPE into the optional arguments

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
NrunsKM = 20;
NiterKM = 500;

% How to initialize GMM
initializationMode = 'kmean';
%initializationMode = 'random';

NrunsGM = 1;
NiterGM = 1000;
check_covar = 0;

% Debug figures
debug = 0;

% Covariance form (netlab vocabulary):
covartype = 'full';

%- Load user options:
if nargin-2 > 0
	if mod(nargin-2,2) ~= 0
		error('Options must come in pairs: OPTION,VALUE')
	end% if 
	for in = 1 : 2 : nargin-2
		eval(sprintf('%s = varargin{in+1};',varargin{in}));
	end% for in	
	clear in
end% if

%- Variables
[dimension np] = size(matrix);

%- Initialize the Gaussian mixture model 
mix0 = gmm(dimension, number_classes, covartype); % init a model structure

switch initializationMode
	case 'kmean' %-- Initialise GMM with a Kmean
		% Update mix0 centres and covariances with KMean results using the gmminit routine
		
		% We run it 'NrunsKM' times and select the most likely model
		% defined as the one with the minimum error.
		% Error is the total squared distance from cluster centres		
		options = foptions;
		options(1)  = 0; % Debug display
		options(2)  = 1e-6; % absolute precision required for the value of CENTRES at the solution. 
							% If the absolute difference between the values of CENTRES between two 
							% successive steps is less than OPTIONS(2), then this condition is satisfied.
		options(3)  = 1e-6; % Difference in successive errors that will stop the iterations (convergence level)
		options(14) = NiterKM;	% Max number of iterations of k-means
		for ir = 1 : NrunsKM
			[mixKM(ir), errlogKM(:,ir), opt, cvgcause{ir}] = gmminit(mix0, matrix', options);
			err(ir) = opt(8); % total squared distance from cluster centres
		end% for ir
		irbestKM = find(err==min(err));

		if debug
			ffland;  iw=1;jw=2;ipl=0;
	
			ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
			p=plot(errlogKM,'k');
			set(p(irbestKM),'color','r','linewidth',2)
			set(gca,'xlim',[0 NiterKM])
			grid on,box on
			xlabel('iterations of k-means');
			ylabel('error');
			title(sprintf('Error'),'fontweight','bold','fontsize',14)
	
			ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
			v = errlogKM;
			v(v==0) = NaN;	
			v = -diff(v,[],1);
			p = plot(v,'k');
			set(p(irbestKM),'color','r','linewidth',2)
			set(gca,'ylim',[1e-12 max(get(gca,'ylim'))])
			set(gca,'xlim',[0 NiterKM])
			hline(options(3))
			set(gca,'yscale','log')	
			grid on,box on
			xlabel('iterations of k-means');
			ylabel('rate of change in error');
			title(sprintf('Rate of change of error'),'fontweight','bold','fontsize',14)

			tt = suptitle(sprintf('KMeans (N-runs=%i)\nMinimal total squared distance from cluster centres',NrunsKM));
			set(tt,'fontweight','bold','fontsize',14)

			clear iw jw ipl p v
		end% if

		if debug & length(irbestKM)>1
			disp('Found more than one KMean with the same minimal error');
		end% if 

		if length(irbestKM)>1
			irbestKM = irbestKM(1);
		end% if 

		%--- Select the KMean solution with the smallest error:
		irbestKM = irbestKM(1);
		errlogKM = errlogKM(:,irbestKM);
		mix0     = mixKM(irbestKM);
		%clear ir opt cvgcause err irbest mix errlog
		
	case 'random' %-- Initialise GMM randomly with data
		% Update mix0 centres and covariances with random stats of random slice of the dataset
	
		% Init with NrunsGM random structures
		for ir = 1 : NrunsGM
			mix0(ir) = gmm(dimension, number_classes, covartype); 
			% Take number_classes random samples within the dataset:
			p = randperm(np);
			p = p(1:number_classes*fix(np./number_classes));
			p = reshape(p,[number_classes fix(np./number_classes)]);
			% Update mix0(ir) centres and covariances:
			centres = mix0(ir).centres;
			covars  = mix0(ir).covars;
			for ik = 1 : number_classes
				slice = matrix(:,p(ik,:));
				
				m = mean(slice,2);
				centres(ik,:) = m;
				
				s = var(slice,[],2);
				switch covartype
					case 'spherical'
						covars(1,ik) = mean(s);
					case 'diag'
						covars(ik,:) = s;
					case 'full'
						covars(:,:,ik) = diag(s);
				end% switch 
				
			end% for ik
			mix0(ir).centres = centres;
			mix0(ir).covars  = covars;
		end% for ir
		errlogKM = NaN;
		
		% Scala spark random init method:
		  % case None => {
		  % 	        val samples = breezeData.takeSample(withReplacement = true, k * nSamples, seed)
		  % 	        (
		  %       Array.fill(k)(1.0 / k), 
		  % 				Array.tabulate(k) 
		  % 				  { 
		  % 					i => 
		  % 						val slice = samples.view(i * nSamples, (i + 1) * nSamples)
		  %         		new MultivariateGaussian(vectorMean(slice), initCovariance(slice))
		  %         }
		  % 	        )
		
		
end% switch 

%- Train GMM with EM

%-- Options
options = foptions;
options(1)  = 0;		% Prints out error values.
%options(3)  = 0;		% Test log likelihood for termination if options(3) is not null
options(3)  = 1e-6;		% Test log likelihood for termination (llh rate of change)
options(5)  = check_covar; % Ensure that covariances don't collapse
options(14) = NiterGM;		% Number of iterations.
options(19) = 0; % Debug by plot in live EM iterations

%-- Training
for ir = 1 : length(mix0)
	[mixGM(ir), opt, errlogGM(:,ir)] = gmmem(mix0(ir), matrix', options, [], []);
	llh(ir) = opt(8); % -sum(log(gmmprob(mix, x)));
end% for ir
irbestGM = find(llh==min(llh));
if debug & length(irbestGM)>1
	disp('Found more than one GMM with the same minimal error');
end% if

%-- Debug figure for GMM
if debug
	ffland;  iw=2;jw=1;ipl=0;
	
	ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
	p = plot(errlogGM,'k-');
	set(p(irbestGM),'color','r','linewidth',2)
	set(gca,'xlim',[0 NiterGM])
	grid on,box on
	xlabel('iterations of GMM');
	ylabel('error');
	title(sprintf('Error'),'fontweight','bold','fontsize',14)
		
	ipl=ipl+1;subp(ipl)=subplot(iw,jw,ipl);hold on
	v = errlogGM;
	v(v==0) = NaN;	
	v = -diff(v,[],1);
	p = plot(v,'k-');
	set(p(irbestGM),'color','r','linewidth',2)
	set(gca,'ylim',[1e-12 max(get(gca,'ylim'))])
	set(gca,'xlim',[0 NiterGM])
	hline(options(3))
	set(gca,'yscale','log')
	grid on,box on
	xlabel('iterations of GMM');
	ylabel('rate of change in error');
	title(sprintf('Rate of change of error'),'fontweight','bold','fontsize',14)
	
	tt = suptitle(sprintf('GMM\nError value is the negative log likelihood of data'));
	set(tt,'fontweight','bold','fontsize',14)
	
end% if 

%-- Select best models for outputs
irbestGM = irbestGM(1);
mixGM = mixGM(irbestGM);
errlogGM = errlogGM(:,irbestGM);
log_likelihood = llh(irbestGM);

%- Compute other usefull metrics for outputs
post = gmmpost(mixGM,matrix');
prob = gmmprob(mixGM,matrix');
acti = gmmactiv(mixGM,matrix');
[~,main] = max(post,[],2);

end %functionem_gmm_v2