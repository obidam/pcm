# -*coding: UTF-8 -*-
#
# To create the sample dataset:
#   ncroot = '/Users/gmaze/data/ARGO/copoda_db/setup_H/db_thd_config6_last/gmm'
#   ncfile = 'NATL_HOMOGENEOUS_variables_7subset_1.nc'
#   dtrain = xr.open_mfdataset(os.path.join(ncroot, ncfile))
#   dtrain.sel(DEPTH=slice(0,-100),N_PROF=slice(0,100)).to_netcdf('argo_sample_test.nc')
#
__author__ = 'gmaze@ifremer.fr'

import unittest
import xarray as xr
import numpy as np
import sklearn as sk
from sklearn.utils.estimator_checks import check_estimator
import pcmpack as PCM
import sys
sys.dont_write_bytecode = True  # This prevent to create .pyc files, useful for debug and reload

class TestV2(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dummydataset = xr.open_mfdataset('argo_sample_test.nc')
        cls.X = cls.dummydataset.TEMP
        cls.Z = cls.dummydataset.DEPTH
        cls.pcm = PCM.v2(3, cls.dummydataset.DEPTH)

    @classmethod
    def tearDownClass(cls):
        del cls.dummydataset
        del cls.X
        del cls.Z
        del cls.pcm

    def test_fit(self):
        self.pcm = self.pcm.fit(self.X, self.Z)
        self.assertTrue(self.pcm._trained)

    def test_predict(self):
        labels = self.pcm.predict(self.X, self.Z)
        self.assertIsNone(sk.utils.assert_all_finite(labels))
        self.assertFalse(np.any(np.logical_or(labels < 0, labels > self.pcm.K)))

    def test_fit_predict(self):
        labels = self.pcm.fit_predict(self.X, self.Z)
        self.assertIsNone(sk.utils.assert_all_finite(labels))
        self.assertFalse(np.any(np.logical_or(labels < 0, labels > self.pcm.K)))

    def test_predict_proba(self):
        post = self.pcm.predict(self.X, self.Z)
        self.assertTrue(np.any(np.logical_and(post >= 0, post <= 1)))

class TestV3(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dummydataset = xr.open_mfdataset('argo_sample_test.nc')
        cls.X = cls.dummydataset.TEMP
        cls.Z = cls.dummydataset.DEPTH
        cls.pcm = PCM.v3(3, cls.dummydataset.DEPTH)

    @classmethod
    def tearDownClass(cls):
        del cls.dummydataset
        del cls.X
        del cls.Z
        del cls.pcm

    def test_fit(self):
        self.pcm = self.pcm.fit(self.X, self.Z)
        self.assertTrue(self.pcm._trained)

    def test_predict(self):
        labels = self.pcm.predict(self.X, self.Z)
        self.assertIsInstance(labels, xr.core.dataarray.DataArray)
        self.assertIsNone(sk.utils.assert_all_finite(labels))
        self.assertFalse(np.any(np.logical_or(labels < 0, labels > self.pcm.K)))

    def test_fit_predict(self):
        labels = self.pcm.fit_predict(self.X, self.Z)
        self.assertIsInstance(labels, xr.core.dataarray.DataArray)
        self.assertIsNone(sk.utils.assert_all_finite(labels))

    def test_predict_proba(self):
        post = self.pcm.predict(self.X, self.Z)
        self.assertIsInstance(post, xr.core.dataarray.DataArray)

class TestV4(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dummydataset = xr.open_mfdataset('argo_sample_test.nc')
        cls.X = cls.dummydataset.TEMP
        cls.Z = cls.dummydataset.DEPTH
        cls.pcm = PCM.v4(3, cls.dummydataset.DEPTH)

    @classmethod
    def tearDownClass(cls):
        del cls.dummydataset
        del cls.X
        del cls.Z
        del cls.pcm

    def test_scikitestimator(self):
        check_estimator(PCM.v4)

    def test_fit(self):
        self.pcm = self.pcm.fit(self.X, axis=self.Z)
        self.assertTrue(self.pcm._trained)

    def test_predict(self):
        labels = self.pcm.predict(self.X, axis=self.Z)
        self.assertIsInstance(labels, xr.core.dataarray.DataArray)
        self.assertIsNone(sk.utils.assert_all_finite(labels))
        self.assertFalse(np.any(np.logical_or(labels < 0, labels > self.pcm.K)))

    # def test_fit_predict(self):
    #     labels = self.pcm.fit_predict(self.X, axis=self.Z)
    #     self.assertIsInstance(labels, xr.core.dataarray.DataArray)
    #     self.assertIsNone(sk.utils.assert_all_finite(labels))
    #
    # def test_predict_proba(self):
    #     post = self.pcm.predict(self.X, axis=self.Z)
    #     self.assertIsInstance(post, xr.core.dataarray.DataArray)

if __name__ == '__main__':
    unittest.main()