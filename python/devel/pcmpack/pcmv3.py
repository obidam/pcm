#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Properties and methods similar to v2
But methods return xarray rather than arrays

Created on 2017/09/26
@author: gmaze
"""

import numpy as np
import xarray as xr
from scipy import interpolate
import copy

from sklearn.utils import validation

# http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
from sklearn import preprocessing

# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
from sklearn.decomposition import PCA

# http://scikit-learn.org/stable/modules/mixture.html
from sklearn.mixture import GaussianMixture

# Decorators:
def deprec_mess(fct):
    print "Version 3 is in debug, do not use !"
    return fct

#
@deprec_mess
class PCM:
    """
        Common base class for a Profile Classification Model
    """
    def __init__(self, K, DPTmodel, scaling=1, reduction=1, classif='gmm', COVARTYPE='full', maxvar=99.9, verb=False):
        """Create the PCM instance
        """
        if   scaling==0: with_scaler = 'none'; with_mean=False; with_std = False
        elif scaling==1: with_scaler = 'normal';  with_mean=True; with_std = True
        elif scaling==2: with_scaler = 'center';  with_mean=True; with_std = False
        else: raise NameError('scaling must be 0, 1 or 2')
        
        if   reduction==0: with_reducer = False
        elif reduction==1: with_reducer = True
        else: raise NameError('reduction must be 0 or 1')
        
        if classif=='gmm': with_classifier = 'gmm';
        else: raise NameError("classifier must be 'gmm' (no other methods at this time)")
        
        self._props = {'K': np.int(K),
                        'llh': None,
                        'COVARTYPE': COVARTYPE,
                        'with_scaler': with_scaler,
                        'with_reducer': with_reducer,
                        'with_classifier': with_classifier,
                        'maxvar': maxvar,
                        'DPTmodel': np.float32(DPTmodel)}
        self._trained = False #todo _trained is a property, should be set/get with a decorator
        self._verb = verb #todo _verb is a property, should be set/get with a decorator
        
        self._interpoler = self.__Interp(self._props['DPTmodel'])
        
        self._scaler = preprocessing.StandardScaler(with_mean=with_mean,
                                                    with_std=with_std)
        self._reducer = PCA(n_components=self._props['maxvar']/100,
                           svd_solver='full')
        self._classifier = GaussianMixture(n_components=self._props['K'],
                                          covariance_type=self._props['COVARTYPE'],
                                          init_params='kmeans',
                                          max_iter=1000,
                                          tol=1e-6)
        self._version = '0.3'

    def __call__(self, **kwargs):
        self.__init__(**kwargs)
    
    def __iter__(self):
        self.__i = 0
        return self
    
    def next(self):
        if self.__i < self.K:
            i = self.__i
            self.__i += 1
            return i
        else:
            raise StopIteration()

    class __Interp:
        """ Internal machinery for the interpolation of vertical profiles
            
            This class is called once at PCM instance initialisation and
            whenever data to be classified are not on the PCM vertical axis.
            
        """
        def __init__(self,DPTmodel):
            self.zi = DPTmodel
            self.doINTERPz = False
        
        def isnecessary(self,C,z):
            """Check wether or not the input data vertical axis is different
                from the PCM one, if not, avoid interpolation
            """
            #todo We should be smarter and recognize occurences of z in DPTmodel
            # or viceversa in order to limit interpolation as much as possible !
            z = np.float32(z)
            #self.doINTERPz = not np.array_equal(self.zi,z)
            self.doINTERPz = not np.array_equiv(self.zi,z)
            return self.doINTERPz
        
        def mix(self,x):
            """ 
                Homogeneize the upper water column: 
                Set 1st nan value to the first non-NaN value
            """
            #izmixed = np.argwhere(np.isnan(x))
            izok = np.where(~np.isnan(x))[0][0]
            #x[izmixed] = x[izok]
            x[0] = x[izok]
            return x;
        
        def fit_transform(self,C,z):
            """
                Interpolate data on the PCM vertical axis
            """
            if (self.isnecessary(C,z)):
                [Np, Nz] = C.shape            
                # Possibly Create a mixed layer for the interpolation to work 
                # smoothly at the surface
                if ((z[0]<0.) & (self.zi[0] == 0.)):
                    z = np.concatenate((np.zeros(1),z))
                    x = np.empty((Np,1))
                    x.fill(np.nan)
                    C = np.concatenate((x,C),axis=1)
                    np.apply_along_axis(self.mix,1,C)
                # Linear interpolation of profiles onto the model grid:
                #f = interpolate.interp2d(z, np.arange(Np), C, kind='cubic')
                f = interpolate.interp2d(z, np.arange(Np), C, kind='linear')
                C = f(self.zi, np.arange(Np))
            return C

    @property
    def K(self):
        """Return the number of class K in the PCM"""
        return self._props['K']

    def display(self, deep=False):
        """Display detailled parameters of the PCM
            This is not a get_params because it doesn't return a dictionnary
            Set Boolean option 'deep' to True for all properties display
        """
        summary = [("<pcm '%s' (K: %i, Z: %i)>")%(self._props['with_classifier'],self._props['K'],self._props['DPTmodel'].size)]
        
        # PCM core properties:
        prop_info = ('Number of class: %i') % self._props['K']
        summary.append(prop_info)
        
        # prop_info = ('Vertical axis: %s') % self._props['DPTmodel']
        prop_info = ('Vertical axis: [%s, ..., %s]') % (repr(self._props['DPTmodel'][0]),repr(self._props['DPTmodel'][-1]))
        summary.append(prop_info)
        
        prop_info = ('Trained: %r') % self._trained
        summary.append(prop_info)
        
        # PCM workflow parameters:
        prop_info = ('Vertical Interpolation: %r') % self._interpoler.doINTERPz
        summary.append(prop_info)    
        summary.append("\t Interpoler: %s"%(type(self._interpoler)))
        
        prop_info = ('Sample Scaling: %r') % self._props['with_scaler']
        summary.append(prop_info)
        summary.append("\t Scaler: %s"%(type(self._scaler)))
        
        if (deep):
            summary.append("\t Scaler properties:")
            d = self._scaler.get_params(deep=deep)
            for p in d: summary.append(("\t\t %s: %r")%(p,d[p]))
        
        prop_info = ('Dimensionality Reduction: %r') % self._props['with_reducer']
        summary.append(prop_info)       
        summary.append("\t Reducer: %s"%(type(self._reducer)))
        #prop_info = ('\t Maximum Variance: %0.2f%%') % self._props['maxvar']
        #summary.append(prop_info) 
        
        if (deep):
            summary.append("\t Reducer properties:")
            d = self._reducer.get_params(deep=deep)
            for p in d: summary.append(("\t\t %s: %r")%(p,d[p]))
        
        prop_info = ('Classification: %r') % self._props['with_classifier']
        summary.append(prop_info) 
        summary.append("\t Classifier: %s"%(type(self._classifier)))
        #prop_info = ('GMM covariance type: %s') % self._props['COVARTYPE']
        #summary.append(prop_info)
        if (self._trained):
            prop_info = ('\t log likelihood: %f') % self._props['llh']
            summary.append(prop_info)
        
        if (deep):
            summary.append("\t Classifier properties:")
            d = self._classifier.get_params(deep=deep)
            for p in d: summary.append(("\t\t %s: %r")%(p,d[p]))
        
        # Done
        return '\n'.join(summary)
    
    def __repr__(self):
        return self.display(deep=self._verb)
    
    def copy(self):
        """Return a deep copy of the PCM instance"""
        return copy.deepcopy(self)

    def preprocessing(self, X, Z):
        """"Pre-process data for classification

            Preprocessing steps:
                interpolation,
                scaling,
                reduction.

            Parameters
            ----------
            X : array-like, shape (N_p=n_samples, N_z=n_features)
                List of N_z-dimensional data profile. Each row
                corresponds to a single profile.
            Z: array-like, shape (N_z=n_features,)
                Vertical axis of profiles

            Returns
            -------
            X : array-like, shape (N_p=n_samples, n_reduced_scaled_interpolated_features)
                List of profiles pre-processed for classification
        """

        # INTERPOLATION:
        X = self._interpoler.fit_transform(X, Z)

        # SCALING:
        self._scaler.fit(X)
        X = self._scaler.transform(X)

        # REDUCTION:
        if self._props['with_reducer']:
            self._reducer.fit(X)
            X = self._reducer.transform(X)

        # Output:
        return X

    def fit(self, X, Z):
        """Estimate PCM parameters
           
            For a PCM, the fit method consists in the following operations:
                - interpolation to the Depth levels of the model
                - scaling
                - reduction
                - estimate classifier parameters

            Parameters
            ----------
            X : array-like, shape (N_p=n_samples, N_z=n_features)
                List of N_z-dimensional data profile. Each row
                corresponds to a single profile.
            Z: array-like, shape (N_z=n_features,)
                Vertical axis of profiles

            Returns
            -------
            self
        """
        # PRE-PROCESSING:
        X = self.preprocessing(X, Z)
        
        # CLASSIFICATION-MODEL TRAINING:
        self._classifier.fit(X)
        self._props['llh'] = self._classifier.score(X)
        
        # Done:
        self._trained = True
        return self
    
    def predict(self, X, Z):
        """Predict the labels for the profile samples in X using trained PCM

           This method add these properties to the PCM data property:
              llh: The log likelihood of the model with regard to new data

            Parameters
            ----------
            X : array-like, shape (N_p=n_samples, N_z=n_features)
                List of N_z-dimensional data profile. Each row
                corresponds to a single profile.

            Returns
            -------
            labels : xarray.DataArray, shape (N_p = n_samples)
                Component labels.
        """
        # if not self._trained:
        #     raise ValueError("Can't predict before fitting !")
        validation.check_is_fitted(self, '_trained',msg="This %(name)s instance is not fitted yet. Call ‘fit’ with appropriate arguments before using this method.")

        # PRE-PROCESSING:
        X = self.preprocessing(X, Z)

        # CLASSIFICATION PREDICTION:
        labels = self._classifier.predict(X)
        self._props['llh'] = self._classifier.score(X)

        # Prepare xarray for output:
        labels = xr.DataArray(labels, dims='samples', name='LABELS')
        labels.attrs['llh'] = self._props['llh']

        # done:
        return labels
    
    def fit_predict(self, X, Z):
        """Estimate PCM parameters and predict classes
           
            Train a PCM and predict classes in a single step
           
            This method adds two properties to the PCM instance:
                LABELS: The label (class prediction)
                llh: The log likelihood of the model with regard to new data
        """
        # PRE-PROCESSING:
        X = self.preprocessing(X, Z)
        
        # CLASSIFICATION-MODEL TRAINING:
        self._classifier.fit(X)
        
        # CLASSIFICATION PREDICTION:
        labels = self._classifier.predict(X)
        self._props['llh'] = self._classifier.score(X)

        # Prepare xarray for output:
        labels = xr.DataArray(labels, dims='samples', name='LABELS')
        labels.attrs['llh'] = self._props['llh']

        # Done:
        self._trained = True
        return labels
    
    def predict_proba(self, X, Z):
        """Predict posterior probability of each component given the data
           
           This method adds these properties to the PCM instance:
               llh: The log likelihood of the model with regard to new data

            Returns
            -------
            resp : array, shape (n_samples, n_components)
                Returns the probability each Gaussian (state) in
                the model given each sample.
        """
        if (not self._trained):
            raise ValueError("Can't predict before fitting !")
        
        # PRE-PROCESSING:
        X = self.preprocessing(X, Z)
        
        # CLASSIFICATION PREDICTION:
        post = self._classifier.predict_proba(X)
        self._props['llh'] = self._classifier.score(X)

        # Prepare xarray for output:
        post = xr.DataArray(post, dims={'samples','components'}, name='POST')
        post.attrs['llh'] = self._props['llh']

        # done:
        return post

    def quant(self, X, Z=None, labels=None, q=[0.05, 0.5, 0.95], verb=False):
        """Compute the qth quantiles of the data for each PCM component.

            Usage A:
                pcm.quant(X, labels=L, q=[0.05,0.5,0.95])
                    This will use labels L do compute component percentiles of X
            Usage B:
                pcm.quant(X, Z=DEPTH, q=[0.05,0.5,0.95])
                    This will classify data in X and then compute percentiles
                    Be careful, if you re-fit a model, you may not end up with something coherent
                    from previous calculation of labels and posteriors, as components will show
                    up in different orders

            Parameters
            ----------
            X : array-like, shape (N_p=n_samples, N_z=n_features)
                List of N_z-dimensional data profile. Each row
                corresponds to a single profile.
            labels: array, shape (N_p=n_samples,)
                Component labels.
            q: float in the range of [0,1] (or sequence of floats), shape (n_quantiles,1)
                Quantile(s) to compute, which must be between 0 and 1 inclusive.

            Returns
            -------
            Q : xarray.DataArray, shape (K, n_quantiles, N_z=n_features)

        """
        if labels == None:
            labels = self.fit_predict(X,Z)
        elif Z == None:
            if not self._trained:
                raise ValueError("Can't compute quant without a fitted model !")

        #
        if (not isinstance(X,xr.core.dataarray.DataArray)):
            XR = xr.DataArray(X, dims=['samples', 'features'])
        else:
            XR = xr.DataArray(X.values, dims=['samples', 'features'])

        if (not isinstance(labels,xr.core.dataarray.DataArray)):
            LR = xr.DataArray(labels, dims=['samples'])
        else:
            LR = xr.DataArray(labels.values, dims=['samples'])

        DS = xr.Dataset({'DATA': XR, 'LABELS': LR})
        varname = 'DATA'
        Q = [] # list of results
        for label, group in DS.groupby('LABELS'):
            if verb:
                print ("Using %0d profiles of %s in class %i") % (group['samples'].shape[0], varname, label)
            quant = group[varname].quantile(q, dim='samples')
            Q.append(quant)
        Q = xr.concat(Q, dim='components') # Transform the list into a DataArray
        Q.name = 'QUANTILES'

        # Done:
        return Q