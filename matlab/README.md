How to use the Matlab PCM toolbox ?
===================================
Created: 2017-12-06, Guillaume Maze (Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

# Install

## Get the source code
The last version of the toolbox from the git repository:
```
    git clone https://github.com/obidam/pcm.git
```

## Dependencies
At this time, the PCM Matlab toolbox relies on the [Netlab Neural Network Software, version 3.3](http://www.aston.ac.uk/eas/research/groups/ncrg/resources/netlab/). But it is included in the PCM distribution under the `src` folder, so you don't have to install it.

## Install
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

Also note that the toolbox works with vertical axis with negative values and oriented from the surface toward the bottom.

```matlab
% Dummy dataset:
N_SAMPLE = 1000;
N_LEVEL = 50;
temp = sort(rand(N_LEVEL, N_SAMPLE)*10+10, 1, 'descend');
dpt = linspace(0, -1000, N_LEVEL)';

% Plot one random profile with: 
plot(temp(:, randi([1, N_SAMPLE],1)),dpt)
```

Now we can train a PCM on this data, using the same vertical axis:
```matlab	
K = 4; % Number of class to fit
PCM = pcmtrain(temp, K, 'full', dpt, 'maxvar',inf);
% 'full' defines the Gaussian class covariance matrix shape
% 'maxvar'=inf is the maximum variance to be retained during compression step, it is an option and 
%	is used here only because we have random data hard to reduce
```

The trained Profile Classification Model is in the PCM structure where all properties should be quite self-explanatory. 
But no direct use of the PCM property should be done.
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

Now you can classify your dataset:
```matlab
[POST LABEL] = pcmpredict(PCM, dpt, temp)
```
The 'POST' array is the class component probability for each data, 'LABEL' is the class maximising 'POST' for each data.

## Saving/loading a PCM:
The toolbox comes with handy functions to save and load the PCM in a netcfd format.

To save a PCM:
```matlab
	pcmsave('my_first_pcm.nc', PCM);
```

To load PCM:
```matlab
	PCM = pcmload('my_first_pcm.nc');
```

The North-Atlantic PCM trained on Argo data used in Maze et al (Pr.Oc., 2017) is published here:
    http://doi.org/10.17882/47106

## Classify a new profile

It is important to note that the PCM has its own vertical axis, so that any new data defined on another vertical axis can
be predicted with a PCM that was possibility train on another dataset:
```matlab
    % Load a profile classification model structure:
    PCM = pcmload('Argo-NATL-PCM-model.nc'); % Netcdf file can be downloaded here: http://doi.org/10.17882/47106
    % This PCM is working with profiles on 0/1405m

    % Load an Argo profile (using the LOPS-Argo Matlab library):
    [Co,Dim] = libargo.read_netcdf_allthefile('6900828_prof.nc');
    pres = Co.pres_adjusted.data(10,:);
    dpt  = -(pres(:)*1.019716); % % Convert pressure to depth (rough approximation for demonstration purpose only)
    temp = Co.temp_adjusted.data(10,:)';

    % Classify the new Argo profile:
    [POST LABEL] = pcmpredict(PCM, dpt, temp) % Note that dpt and temp MUST be (z,:)
```

This example uses the LOPS-Argo Matlab toolbox to load Argo netcdf data. **The PCM toolbox does not depend on it, this 
is only for examples**. The toolbox is available here:
```
svn checkout https://forge.ifremer.fr/svn/lpoargo/projects/matlab forge-lpoargo
```
and add it to your path:
```matlab
    addpath('forge-lpoargo');
```

## Working with gridded fields 

The PCM matlab toolbox core functions 'train' and 'predict' comes with handy counterparts to work with gridded fields 
rather than collection of profiles.

For the example below, we'll work with one monthly field of the ISAS dataset (you can get it here: http://doi.org/10.17882/45945)
```matlab
% Load ISAS snapshot
ncfile = 'ISAS13_20121215_fld_TEMP.nc';
temp = ncread(ncfile,'TEMP'); % lon/lat/dpt
temp = permute(temp,[3 2 1]); % Because PCM works in dpt/lat/lon
dpth = -double(ncread(ncfile,'depth')); % Because PCM work with negatively defined downward depth axis
lat  = ncread(ncfile,'latitude');
lon  = ncread(ncfile,'longitude');
LSmask = squeeze(double(~isnan(temp(1,:,:)))); % Land=0, Sea=1, mask
```

### Train a PCM on a gridded field

Note that before training the PCM, we need to create a 2D mask defining the location of profiles to be classified (they 
should be plain, without NaNs). This is easily done with the ```createGMMmask``` function.

```matlab
DPTmodel = dpth(find(dpth>=-1000));
MASK = createGMMmask(LSmask, temp, dpth, 'Zmin', min(DPTmodel));
PCM = pcmtrainlatlon(MASK, temp, 8, 'full', DPTmodel, 'DPT', dpth);
```

Note that in the above, we made use of the possibility to train and work with a PCM that is not on the same vertical 
axis of the data.

### Classify a gridded field

Finally, you can classify the gridded field:
```matlab
% Classify the gridded temperature ISAS field:
[POST LABELS ACTI PROB] = pcmpredictlatlon(MASK, PCM, dpth, temp);

% Plot the results:
figure;
subplot(2,1,1);pcolor(lon,lat,LABELS);shading flat
caxis([1 PCM.K]);colormap(jet(PCM.K))

subplot(2,1,2);pcolor(lon,lat,squeeze(POST(:,:,1)));shading flat
caxis([0 1]);
```

In the case where you want to ensure the MASK corresponds to the PCM depth axis, you can update it this way:
```matlab
MASK = createGMMmask(LSmask, temp, dpth, 'Zmin', min(PCM.DPTmodel));
```
