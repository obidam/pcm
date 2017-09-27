#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Try to make pcmv3 a valid scikit-learn estimator

The main motivation to make our class compatible to the scikit-learn estimator interface
is to use it together with model evaluation and selection tools such as:
    model_selection.GridSearchCV
    pipeline.Pipeline

Help:
    http://scikit-learn.org/stable/developers/contributing.html#rolling-your-own-estimator

Created on 2017/09/26
@author: gmaze
"""

import numpy as np
import xarray as xr
from scipy import interpolate
from sklearn.base import BaseEstimator, ClassifierMixin
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
    print "Version 4 is in debug, do not use !"
    return fct

#
@deprec_mess
class PCM(BaseEstimator, ClassifierMixin):
    """Profile Classification Model

        Parameters
        ----------

        Methods
        -------

        Examples
        --------

    """
    def __init__(self, n_components=1, axis=9999, scaling=1, reduction=1, classifier='gmm', COVARTYPE='full', maxvar=99.9, verb=False):
        """Create the PCM instance
        """
        if   scaling==0: with_scaler = 'none'; with_mean=False; with_std = False
        elif scaling==1: with_scaler = 'normal';  with_mean=True; with_std = True
        elif scaling==2: with_scaler = 'center';  with_mean=True; with_std = False
        else: raise NameError('scaling must be 0, 1 or 2')
        
        if   reduction==0: with_reducer = False
        elif reduction==1: with_reducer = True
        else: raise NameError('reduction must be 0 or 1')
        
        if classifier=='gmm': with_classifier = 'gmm';
        else: raise NameError("classifier must be 'gmm' (no other methods at this time)")
        
        self._props = {'K': np.int(n_components),
                        'llh': None,
                        'COVARTYPE': COVARTYPE,
                        'with_scaler': with_scaler,
                        'with_reducer': with_reducer,
                        'with_classifier': with_classifier,
                        'maxvar': maxvar,
                        'DPTmodel': np.float32(axis)}
        self._trained = False #todo _trained is a property, should be set/get with a decorator
        self._verb = verb #todo _verb is a property, should be set/get with a decorator
        self._version = '0.4'

    # def __call__(self, **kwargs):
    #     self.__init__(**kwargs)

    def get_params(self, deep=True):
        # suppose this estimator has parameters "alpha" and "recursive"
        return self._props

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            self._props[parameter] = value
        return self

    def set_config(self, **kargs):
        """Set-up all processing steps according to PCM properties"""

        self.set_params(kargs)

        self._interpoler = self.__Interp(self._props['DPTmodel'])

        self._scaler = preprocessing.StandardScaler(with_mean=with_mean,
                                                    with_std=with_std)
        self._reducer = PCA(n_components=self._props['maxvar'] / 100,
                            svd_solver='full')
        self._classifier = GaussianMixture(n_components=self._props['K'],
                                           covariance_type=self._props['COVARTYPE'],
                                           init_params='kmeans',
                                           max_iter=1000,
                                           tol=1e-6)
        return self

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
    
    # def __repr__(self):
    #     return self.display(deep=self._verb)

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

    def fit(self, X, y=None, axis=None, **kargs):
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
            axis: array-like, shape (N_z=n_features,)
                Vertical axis of profiles

            Returns
            -------
            self
        """
        self.set_config(kargs)

        # PRE-PROCESSING:
        X = self.preprocessing(X, axis)
        
        # CLASSIFICATION-MODEL TRAINING:
        self._classifier.fit(X)
        self._props['llh'] = self._classifier.score(X)
        
        # Done:
        self._trained = True
        return self

    def score(self, X, y=None, axis=None):
        """Compute the per-sample average log-likelihood of the given data X
        """
        if (not self._trained):
            raise ValueError("Can't predict before fitting !")

        # PRE-PROCESSING:
        X = self.preprocessing(X, axis)

        return self._classifier.score(X)

    def predict(self, X, y=None, axis=None):
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
        X = self.preprocessing(X, axis)

        # CLASSIFICATION PREDICTION:
        labels = self._classifier.predict(X)
        self._props['llh'] = self._classifier.score(X)

        # Prepare xarray for output:
        labels = xr.DataArray(labels, dims='samples', name='LABELS')
        labels.attrs['llh'] = self._props['llh']

        # done:
        return labels
