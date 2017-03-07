function varargout = get_fuzzy_classprofiles(EOFs,XR,Xref,mix,LABELS,WEIGHTS,varargin)
% get_fuzzy_classprofiles Compute representative class profiles given a GMM and a dataset (fuzzy PDF)
%
% [MS PRC MM, PDFs] = get_fuzzy_classprofiles(EOFS,XR,XREF,MIX,LABELS,WEIGHTS,[PAR,VAL])
%	Compute representative class profiles from a GMM and a dataset using fuzzy classification
% statistics, i.e. we use weight profiles to compute representative profiles.
%
% Inputs:
% - Parameters used to recompose the original dataset:
%	EOFS [Nz,Nc]: Principal component coefficients, ie eigenvectors
% 	XR   [Nc,Np]: Reduced dataset, ie eigenvalues
% 	XREF [Nz,1]:  Reference profile as output by gmmargo_1param_reduce.m
% - Parameters used to identify classes and compute class statistics:
% 	MIX (struct):   A GMModel structure as output by em_gmm.m
% 	LABELS  [Np,1]: GMMmodel labels as output by em_gmm.m. Labels are used to identify
%	                profiles within a class.
% 	WEIGHTS [Np,K]: Define how to weight profiles to compute class statistics. For
%	                instance, it can be GMMmodel posterior values as output by em_gmm.m.
% - Additional parameters [PAR,VAL]:
% 	oLABELS [Np,1]: Effective labels used in plots, as output by rename_labels.m
%	lcomp [integers]: The list of EOFs to use (any integer between 1 and Nc)
% 	pdfs  [bool]: A boolean value to indicate if we should return the PDFs 
%	rg % Range of the histogram x-axis
% 
% Outputs:
%	MS [2,K,Nz]: Empirical class moments (mean/std) from weighted statistics of profiles.
% 	             MS [1,:,:] is the mean
% 	             MS [2,:,:] is the std
% 	             MS [2,:,:] is the skewnewss
%	PRC[3,K,Nz]: Empirical class percentiles (5/50/95) from profiles fuzzy classified.
% 	             PRC [1,:,:] is the 50% percentile (ie median)
% 	             PRC [2,:,:] is the  5% percentile
% 	             PRC [3,:,:] is the 95% percentile
%	MM [2,K,Nz]: Analytical class moments (mean/std) from GMM class centers converted into original data space.
% 	             MM [1,:,:] is the mean
% 	             MM [2,:,:] is the std [TO BE DONE !, return NaN as of now]
% 	PDFs: A structure with the following probability density functions:
% 		PDFS.sample [Nz,n]:   The observed PDF at each level (computed using the complete sample of profiles)
% 		PDFS.sclass [Nz,K,n]: The empirical model PDF at each level for each class (computed using fuzzy classified profiles)
% 		PDFS.recomp [Nz,n]:   The empirical model PDF at each level from the sum of priors weighted pdfs of each class
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-03-31 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

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
n  = 200; % x-axis resolution of the histogram
xl = [-1 1]*4; % Range of the histogram x-axis
rg = linspace(xl(1),xl(2),n);
pdfs = false;

%- Load user parameters:
if nargin > 6
    if mod(nargin,2) ~= 0
        error('Parameters must come in pairs: PAR,VAL');
    end% if
    for in = 1 : 2 : nargin-6
        eval(sprintf('%s = varargin{in+1};',varargin{in}));
    end% for in
    clear in
end% if

%- Update default parameters with user data
Nz = size(EOFs,1); % = length(DPT)
Nc = size(XR,1);
K  = size(WEIGHTS,2);
Np = size(XR,2);

if exist('Xscl_std') % We're going to compute PDF of de-normalized anomalies
	rgscl = rg*xtrm(Xscl_std);
end% if 

% Nb of EOFs to use to recompose observed data from compressed data
if ~exist('lcomp','var')	
	lcomp = 1:size(XR,1); 
else
	if ~isempty(find(lcomp>Nc))
		error('lcomp cannot be larger than the number of EOFs !')
	end% if 
	if length(unique(lcomp))~=length(lcomp)
		error('You cannot use an EOFs more than once to recompose the signal !')
	end% if 
end% if 

if ~exist('oLABELS','var')
	oLABELS = LABELS;
end% if 

%- Decompose the pdf for each class and compute profiles
if pdfs
	PDFd = zeros(Nz,length(rg));   % Directly from the sample
	PDFc = zeros(Nz,length(rg));   % From the sample, but recomposed from the pdf of each class
	PDFk = zeros(Nz,K,length(rg)); % From the sample, details for each class
	dr = diff(rg(1:2));
	if exist('Xscl_std')	
		PDFc_ano = zeros(Nz,length(rgscl));   % From the sample, details for each class
		PDFk_ano = zeros(Nz,K,length(rgscl)); % From the sample, details for each class
	end% if 
end% if 
Cmean_samp = zeros(K,Nz);
Cmedi_samp = zeros(K,Nz);
C05_samp = zeros(K,Nz);
C95_samp = zeros(K,Nz);

for iz = 1 : Nz
	% Expand dataset (if reduced with PCA method):
	ds = EOFs(iz,lcomp)*XR(lcomp,:) + Xref(iz); % Recompose data for this level from EOFs
	
	%- Stats for the entire sample:
	if pdfs
		%-- PDFs
		c = histc(ds,rg); % Count data
		c = c/sum(c*dr);  % Transform histogram into pdf (integral is one, not nb of points)
		PDFd(iz,:) = c;
		%-- mean/std/skewness and 5%/50%/95% percentile profiles
		Smean_samp(iz) = mean(ds);
		Sstd_samp(iz)  = std(ds);
		S05_samp(iz)   = prctile(ds, 5,2);
		Smedi_samp(iz) = prctile(ds,50,2);
		S95_samp(iz)   = prctile(ds,95,2);
	end% if 
	
	%- Stats for each class:
	for ik = 1 : K
		ikmix = unique(LABELS(find(oLABELS==ik)));
		
		% Should we de-norm ? and de-scale ?
		if pdfs
			c = histcw(ds,WEIGHTS(:,ikmix),rg); % Count data weighted by WEIGHTS
			c = c/sum(c*dr);  % Transform weighted histogram into pdf (integral is one, not nb of points)
			PDFk(iz,ik,:) = c;
			PDFc(iz,:) = PDFc(iz,:) + squeeze(mix.priors(ikmix)*PDFk(iz,ik,:))';
			if exist('Xscl_std')
				ds = ds.*Xscl_std(iz);				
				c = histcw(ds,WEIGHTS(:,ikmix),rgscl); % Count data weighted by WEIGHTS				
				c = c/sum(c*diff(rgscl(1:2)));  % Transform histogram into pdf (integral is one, not nb of points)
				PDFk_ano(iz,ik,:) = c;
				PDFc_ano(iz,:) = PDFc_ano(iz,:) + squeeze(mix.priors(ikmix)*PDFk_ano(iz,ik,:))';
			end% if
		end% if 
		
		%-- Compute weighted mean/std profiles for this class:
		Cmean_samp(ik,iz)  = wmean(WEIGHTS(:,ikmix),ds);
		Cstd_samp(ik,iz)   = wstd(WEIGHTS(:,ikmix),ds'); 
		Cskew_samp(ik,iz)  = wskewness(WEIGHTS(:,ikmix),ds'); 
		
		%-- Compute 5%/50%/95% percentile of profiles for this class:
		C05_samp(ik,iz)   = wprctile(ds, 5,WEIGHTS(:,ikmix));
		Cmedi_samp(ik,iz) = wprctile(ds,50,WEIGHTS(:,ikmix));
		C95_samp(ik,iz)   = wprctile(ds,95,WEIGHTS(:,ikmix));
		
	end% for ik
end% for iz

%- Project class model centers on the original dimensional space:
Cmean_model = zeros(K,Nz);
Cstd_model  = zeros(K,Nz);
for ik = 1 : K
	% Identify gmm matrix k index
	ikmix  = unique(LABELS(find(oLABELS==ik))); 

	% Get center in reduced-dimension space:
	cxr = mix.centres(ikmix,:)'; 
	Cmean_model(ik,:) = EOFs(:,lcomp)*cxr(lcomp,:) + Xref;  

	% Get std in reduced-dimension space:
	% TODO: I don't know how to do this yet ! ie I fill with NaNs
	switch mix.covar_type
		case 'full',         cxd = sqrt(diag(mix.covars(:,:,ikmix)));
		case {'diag','mix'}, cxd = Inf;
		case 'spherical',    cxd = Inf;
	end% switch	
	% We reconstruct only using the diagonal of the covariance matrix
%	Cstd_model(ik,:) = (EOFs(:,lcomp)*(cxr(lcomp,:)+cxd(lcomp,:)) - EOFs(:,lcomp)*(cxr(lcomp,:)-cxd(lcomp,:)))/2;

%	Cstd_model(ik,:) = EOFs(:,lcomp)*cxr(lcomp,:) + Xref;  % NOT WORKING !
	Cstd_model(ik,:) = NaN;
end% for ik

%- Outputs
MOMsamp(1,:,:) = Cmean_samp;
MOMsamp(2,:,:) = Cstd_samp;
MOMsamp(3,:,:) = Cskew_samp;

PRCsamp(1,:,:) = Cmedi_samp;
PRCsamp(2,:,:) = C05_samp;
PRCsamp(3,:,:) = C95_samp;

MOMmodel(1,:,:) = Cmean_model;
MOMmodel(2,:,:) = Cstd_model;

varargout(1) = {MOMsamp};
varargout(2) = {PRCsamp};
varargout(3) = {MOMmodel};

if pdfs
% 		PDFS.sample [Nz,n]:   The PDF at each level computed using the complete sample
% 		PDFS.sclass [Nz,K,n]: The PDF at each level for each class, computed using the weighted sample
% 		PDFS.recomp [Nz,n]:   The PDF at each level from the sum of priors weighted pdfs of each class
	PDFS.sample = PDFd;
	PDFS.sclass = PDFk;
	PDFS.recomp = PDFc;
	PDFS.sample_mean = Smean_samp;
	PDFS.sample_std  = Sstd_samp;
	PDFS.sample_perc05 = S05_samp;
	PDFS.sample_perc50 = Smedi_samp;
	PDFS.sample_perc95 = S95_samp;
	if exist('Xscl_std')
		PDFS.sclass_ano = PDFk_ano;
		PDFS.recomp_ano = PDFc_ano;
		PDFS.rgscl = rgscl;
	end% if 
	PDFS.readme = 'All profiles were used and histograms computed using weighted counts';
	varargout(4) = {PDFS};
end% if 

end %functionget_fuzzy_classprofiles