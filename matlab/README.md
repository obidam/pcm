How to use the Matlab PCM toolbox ?
===================================
Created: 2017-12-06, Guillaume Maze (Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

# Download
You can get the last version of the toolbox from the git repository:
```
    git clone https://github.com/obidam/pcm.git
```
# Dependencies
At this time, the PCM Matlab toolbox relies on the [Netlab Neural Network Software, version 3.3](http://www.aston.ac.uk/eas/research/groups/ncrg/resources/netlab/). But it is included in the PCM distribution under the `src` folder, so you don't have to install it.

# Install
Simply ensure that you have the PCM Matlab toolbox (`src`) and its dependencies (`lib`) folders on your Matlab path !

If you used the above checkout, yo have to add to your path:
```matlab
    addpath(fullfile('pcm','matlab','src'));
    addpath(fullfile('pcm','matlab','lib','netlab3_3'));
```

# Usage
Here are scripts to introduce you with the toolbox

## Train a PCM on your data

Let's create a dummy dataset: an array 'temp' of dimension [N_LEVEL, N_SAMPLE], and its vertical axis 'dpt'.

This could be any collection of vertical profiles, as long as it is a plain matrix, without NaNs.

```matlab

	N_SAMPLE = 1000;
	N_LEVEL = 50;
	temp = sort(rand(N_LEVEL, N_SAMPLE)*10+10, 1, 'descend');
	dpt = linspace(0, -1000, N_LEVEL)';
	
	% Plot one random profile with: 
	plot(temp(:, randi([1, N_SAMPLE],1)),dpt)
```

Now we train a PCM on this data, using the same vertical axis:
```matlab	
	K = 4; % Number of class
	[PCM DATAi DATAn DATAr] = pcmtrain(temp, K, 'full', dpt, 'maxvar',inf);
	% 'full' defines the Gaussian class covariance matrix shape
	% 'maxvar'=inf is the maximum variance to be retained during compression step, it is an option and 
	%	is used here only because we have random data hard to reduce
```

The trained Profile Classification Model is in the PCM structure:
```matlab
	>> PCM
	PCM =
	  struct with fields:
	         DPTmodel: [501 double]
	             EOFs: [5050 double]
	                K: 4
	              LLH: 13634.1879135846
	               Np: 1000
	                V: [150 double]
	            X_ave: [501 double]
	            X_ref: [501 double]
	            X_std: [501 double]
	        covarTYPE: 'full'
	         doREDUCE: 1
	           maxvar: 100
	              mix: [11 struct]
	    normalization: 1
	           readme: 'This PCM was created using pcmtrain.m'
```

## Classify a single profile

    % Load a classification model structure:
    PCM = pcmload('Argo-NATL-PCM-model.nc'); % Netcdf file can be downloaded here: http://doi.org/10.17882/47106

    % Load an Argo profile (using the LOPS-Argo Matlab library):
    [Co,Dim] = libargo.read_netcdf_allthefile('6900828_prof.nc');
    pres = Co.pres_adjusted.data(10,:);
    dpt  = -(pres(:)*1.019716); % % Convert pressure to depth (rough approximation for demonstration purpose only)
    temp = Co.temp_adjusted.data(10,:)';

    % Classify the new Argo profile:
    [POST LABEL] = pcmpredict(PCM,dpt,temp) % Note that dpt and temp MUST be (z,:)

Note that the above script use the LOPS-Argo Matlab toolbox to load the Argo netcdf data:

    addpath('forge-lpoargo'); % svn checkout https://forge.ifremer.fr/svn/lpoargo/projects/matlab forge-lpoargo

## Classify a gridded field

    % Load a classification model structure:
    PCM = pcmload('Argo-NATL-PCM-model.nc'); % Netcdf file can be downloaded here: http://doi.org/10.17882/47106

    % Load ISAS snapshot
    ncfile = 'NRTOAGL01_20161015_fld_TEMP.nc'; % ISAS OI DATA can be downloaded here: http://doi.org/10.17882/45945
    temp = ncread(ncfile,'TEMP'); % lon/lat/dpt
    temp = permute(temp,[3 2 1]); % Because PCM works in dpt/lat/lon
    dpth = -double(ncread(ncfile,'depth')); % Because PCM work with negatively defined depth axis
    lat  = ncread(ncfile,'latitude');
    lon = ncread(ncfile,'longitude');

    % Mask data that can be classified (defined over a depth range similar than the PCM):
    MASK = repmat(dpth,[1 length(lat) length(lon)]);
    MASK(isnan(temp)) = 0; MASK = squeeze(min(double(MASK)));
    MASK = double(MASK<min(PCM.DPTmodel));

    % Classify the gridded temperature ISAS field:
    [POST LABELS ACTI PROB] = pcmpredictlatlon(MASK,PCM,dpth,temp);

    % Plot the results:
    figure;
    subplot(2,1,1);pcolor(lon,lat,LABELS);shading flat
    caxis([1 PCM.K]);colormap(jet(PCM.K))

    subplot(2,1,2);pcolor(lon,lat,squeeze(POST(:,:,1)));shading flat
    caxis([0 1]);