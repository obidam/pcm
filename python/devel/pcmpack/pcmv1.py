#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 15:30:02 2017

The goal is to create a class that behaves like any scikit method and can
be used into a Pipeline
(So it needs to return "self")

The structure is taken from the list of methods used in the Matlab implementation:
    - pcmtrain
    - pcmpredict
    
The python class is an instance of the PCM Matlab structure returned by a 
pcmsave/pcmload    

# Set-up and train the classifier:
gmm = GaussianMixture(n_components=K,
                      covariance_type='full',
                      init_params='kmeans',
                      max_iter=1000,
                      tol=1e-6)
gmm.fit(Xr) # Training on reduced data

# Extract GMM parameters:
priors = gmm.weights_ # [K,1]
centers= gmm.means_   # [K,Nc]
covars = gmm.covariances_ # [K,Nc,Nc] if 'full'

# Classify the dataset:
LABELS = gmm.predict(Xr) # [Np,1]
POST   = gmm.predict_proba(Xr) # [Np,Nc]


@author: gmaze
"""

# Import all necessary packages:
import numpy as np

# http://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
from sklearn import preprocessing

# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
from sklearn.decomposition import PCA

# http://scikit-learn.org/stable/modules/mixture.html
from sklearn.mixture import GaussianMixture

from scipy import interpolate

def deprec_mess(fct):
    print "Version 1 is deprecated !"
    return fct

#
@deprec_mess
class PCM:
    """
        Common base class for a Profile Classification Model
    """
        
    def __init__(self, K, DPTmodel, scaling=1, reduction=1, classif='gmm', COVARTYPE='full', maxvar=99.9, verb=False):
        #todo: check inputs validity
        if   (scaling==0): with_scaler = False; with_mean=False; with_std = False;
        elif (scaling==1): with_scaler = True;  with_mean=True; with_std = True;
        elif (scaling==2): with_scaler = True;  with_mean=True; with_std = False;
        else: raise NameError('scaling must be 0, 1 or 2');
        
        if   (reduction==0): with_reducer = False; 
        elif (reduction==1): with_reducer = True;  
        else: raise NameError('reduction must be 0, 1');
        
        if   (classif=='gmm'): with_classifier = 'gmm'; 
        else: raise NameError("classifier must be 'gmm' (no other methods at this time)");
                            
        self._props = {'K':np.int(K),
                      'COVARTYPE':COVARTYPE,
                      'with_scaler': with_scaler,
                      'with_reducer': with_reducer,
                      'with_classifier': with_classifier,
                      'maxvar':maxvar,
                      'DPTmodel':np.float32(DPTmodel)};
        self._trained = False;
        self._verb = verb;
        self.K = self._props['K']
        self._scaler  = preprocessing.StandardScaler(with_mean=with_mean,
                                                    with_std=with_std);
        self._reducer = PCA(n_components=self._props['maxvar']/100,
                           svd_solver='full')
        self._classifier = GaussianMixture(n_components=self._props['K'],
                                          covariance_type=self._props['COVARTYPE'],
                                          init_params='kmeans',
                                          max_iter=1000,
                                          tol=1e-6)
        self.interpoler = PCM.Interp(self._props['DPTmodel'])
        self._version = '0.1'

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

    class Interp:
        """ 
            Internal machinery for the interpolation of vertical profiles
            
            This class is called once at PCM instance initialisation and
            whenever data to be classified are not on the PCM vertical axis.
            
        """
        def __init__(self,DPTmodel):
            self.zi = DPTmodel
            self.doINTERPz = False
            
        def isnecessary(self,C,z):
            """
                Check whether or not the data vertical axis is different
                from the PCM one, if not, avoid interpolation
            """
            z = np.float32(z)
            # self.doINTERPz = not np.array_equal(self.zi,z)
            self.doINTERPz = not np.array_equiv(self.zi,z)
            return self.doINTERPz
            
        def mix(self,x):
            """ 
                Homogenize the upper water column:
                Set 1st nan value to the first non-NaN value
            """
            # izmixed = np.argwhere(np.isnan(x))
            izok = np.where(~np.isnan(x))[0][0]
            # x[izmixed] = x[izok]
            x[0] = x[izok]
            return x
        
        def fit_transform(self,C,z):
            """
                Interpolate data on the PCM vertical axis
            """
            if self.isnecessary(C,z):
                [Np, Nz] = C.shape            
                # Possibly Create a mixed layer for the interpolation to work 
                # smoothly at the surface
                if (z[0]<0.) & (self.zi[0] == 0.):
                    z = np.concatenate((np.zeros(1),z))
                    x = np.empty((Np,1)); x.fill(np.nan)
                    C = np.concatenate((x,C),axis=1)                
                    np.apply_along_axis(self.mix,1,C)
                # Linear interpolation of profiles onto the model grid:
                # f = interpolate.interp2d(z, np.arange(Np), C, kind='cubic')
                f = interpolate.interp2d(z, np.arange(Np), C, kind='linear')
                C = f(self.zi, np.arange(Np))
            return C        
          
    def display(self,deep=True):
        """
            Display detailed parameters of the PCM
            This is not get_params because it should return a dictionnary
        """
        summary = [("<pcm '%s' (K: %i, Z: %i)>")%(self._props['with_classifier'],self._props['K'],self._props['DPTmodel'].size)]
        
        # PCM core properties:
        prop_info = 'Number of class: %i'.format(self._props['K'])
        summary.append(prop_info)
        
        prop_info = ('Vertical axis: %s') % repr(self._props['DPTmodel'])
        summary.append(prop_info)
        
        prop_info = ('Trained: %r') % self._trained        
        summary.append(prop_info)
        
        # PCM workflow parameters:
        prop_info = ('Vertical Interpolation: %r') % self.interpoler.doINTERPz
        summary.append(prop_info)    
        summary.append("\t Interpoler: %s"%(type(self.interpoler)))
    
        prop_info = ('Sample Normalisation: %r') % self._props['with_scaler']
        summary.append(prop_info)
        summary.append("\t Normaliser: %s"%(type(self._scaler)))
    
        if (deep):
            summary.append("\t Normaliser properties:")
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
            prop_info = ('\t log likelihood: %f') % self.llh
            summary.append(prop_info)

        if (deep):            
            summary.append("\t Classifier properties:")
            d = self._classifier.get_params(deep=deep)
            for p in d: summary.append(("\t\t %s: %r")%(p,d[p]))
        
        # Done
        return '\n'.join(summary)
    
    def __repr__(self):
        return self.display(deep=self._verb)
    
    def fit(self,X,Z):
        """
            For a PCM, the fit method consists in the following operations:
                - interpolation to the Depth levels of the model
                - scaling
                - reduction
                - estimate GMM parameters                
        """        
        # CHECK INPUTS
        #todo we should check for errors/inconsistencies in inpts                 
        
        # INTERPOLATION:
        X = self.interpoler.fit_transform(X,Z)

        # SCALING:
        self._scaler.fit(X)
        X = self._scaler.transform(X)

        # REDUCTION:  
        if (self._props['with_reducer']):
            self._reducer.fit(X)
            X = self._reducer.transform(X)

        # CLASSIFICATION-MODEL TRAINING:
        self._classifier.fit(X)
        self.llh = self._classifier.score(X)
        
        # Done:  
        self._trained = True
        return self
    
    def predict(self,X,Z):
        """
            Using the self PCM properties, predict the class of new data
        """
        # CHECK INPUTS
        #todo we should check for errors/inconsistencies in inpts
        print self._trained
        if not self._trained:
            raise NameError("Can't predict before fitting !")

        # INTERPOLATION:
        X = self.interpoler.fit_transform(X,Z)

        # SCALING:
        X = self._scaler.transform(X)

        # REDUCTION:  
        if (self._props['with_reducer']):
            X = self._reducer.transform(X)      

        # CLASSIFICATION PREDICTION:
        self.LABELS = self._classifier.predict(X)
        self.llh = self._classifier.score(X)

        # done:
        return self

    def fit_predict(self,X,Z):
        """
            Train a PCM and predict classes in a single step
        """
        # CHECK INPUTS
        #todo we should check for errors/inconsistencies in inpts                 
        
        # INTERPOLATION:
        X = self.interpoler.fit_transform(X,Z)

        # SCALING:
        self._scaler.fit(X)
        X = self._scaler.transform(X)

        # REDUCTION:     
        if (self._props['with_reducer']):
            self._reducer.fit(X)
            X = self._reducer.transform(X)   

        # CLASSIFICATION-MODEL TRAINING:
        self._classifier.fit(X)
        
        # CLASSIFICATION PREDICTION:
        self.LABELS = self._classifier.predict(X)
        self.llh = self._classifier.score(X)

        # Done:
        self._trained = True
        return self

