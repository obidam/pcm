% Version alpha 0.1 01-03-2017
% 
% 		PCM Toolbox
% 
% Last update: 2017 January 18, 16:31
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017,  OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% 
% 	classpermute                             - Permute class order in a Profile Classification Model (PCM)
% 	covarrotate                              - Rotate a squared covariance matrix
% 	createGMMmask                            - Create a mask from land/sea and 3D field data constraints
% 	demo_gmm_v1                              - Simplest demonstrate of GMM in 2D
% 	demo_gmm_v2                              - Demonstrate GMM in 2D with a T/S diagram
% 	demo_gmm_v3                              - Same as demo_gmm_v2 but with interactive nb of cluster
% 	demo_gmm_v4                              - Another Demonstrate of GMM
% 	em_gmm_v1                                - (deprec) Train a GMM model with EM algorithm and KMean initialisation.
% 	em_gmm_v2                                - (active) Train a GMM model with EM algorithm and KMean initialization.
% 	get_fuzzy_classprofiles                  - Compute representative class profiles given a GMM and a dataset (fuzzy PDF)
% 	get_hard_classprofiles                   - Compute representative class profiles given a GMM and a dataset (hard PDF)
% 	gmmbic                                   - Return the BIC from a GMM and a dataset
% 	pcmdist                                  - Compute the Mahalanobis distance of a collection of profiles to PCM Gaussian
% 	pcmload                                  - Load a Profile Classification Model (PCM) from a netcdf file
% 	pcmoptimk                                - Determine the optimum number of classes K (using BIC)
% 	pcmpdfz                                  - Compute the PDF at each depth levels (PDFz) from a PCM and a given collection of profiles
% 	pcmpdfzlatlon                            - Compute the PDF at each depth levels (PDFz) from a PCM and a given collection of gridded profiles
% 	pcmprctile                               - Compute the percentile at each depth levels from a PCM and a given collection of profiles
% 	pcmprctilelatlon                         - Compute the percentile at each depth levels from a PCM and a given collection of gridded profiles
% 	pcmpredict                               - Classify a collection of profiles with a Profile Classification Model (PCM) based on GMM
% 	pcmpredictlatlon                         - Classify a lat/lon gridded field with Profile Classification Model based on GMM
% 	pcmreduce                                - Reduce a data set using a given PCM
% 	pcmsave                                  - Save a Profile Classification Model (PCM) into a netcdf file
% 	pcmtrain                                 - Train a Profile Classification Model (PCM) based on GMM onto a collection of profiles
% 	pcmtrainlatlon                           - Train a Profile Classification Model (PCM) based on GMM onto a depth/lat/lon gridded field
% 	pcollecnorm                              - Normalize a collection of profiles
% 	plot_GMMellipse                          - Plot Gaussian mode ellipse for a given GMM
% 	reduce_dimensions                        - Reduce the dimension of a dataset with PCA/SVD
% 	rename_labels                            - reorder_labels Re-name labels IDs
% 
%
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