"""
"""
__all__ = ['SpatialLocation']

import numpy as np
import pandas as pd
import healpy as hp
from astropy.table import Table, Column
from astropy.io import fits
from .times import get_jd_from_mjd


class SpatialLocation():
    """

    Parameters
    ----------
    merger_times : sequence of `np.floats`
        times of merger of events
    spatial_fits : sequence of strings 
    """

    def __init__(self, spatial_fits_files, merger_times=None,
                 probability_fractions=None,
                 event_names=None):
        """

        """
        self.fnames = np.ravel(spatial_fits_files)
        if merger_times is None:
            merger_times = self.get_merger_times(self.fnames)
        
        self.times = np.ravel(merger_times)

        if event_names is None:
            self.event_names = np.arange(len(merger_times))
        else:
            self.event_names = np.ravel(event_names)



        if probability_fractions is None:
            probability_fractions = np.ones_like(self.times) * 0.9

        self.prob_fractions_threshold = np.ravel(probability_fractions)
        self._sky_probs = None
        self._skypix = None


    @staticmethod
    def get_merger_times(fnames):
        """ return a list of merger times in Julian Date from the
        LIGO files

        Parameters
        ---------
        fnames : sequence of file names of LIGO event inferred from the
            data

        Returns
        -------
        sequence of julian days of merger times.
        """

        jds = list()
        for fname in fnames:
            data = fits.open(fname)
            mjd = data[1].header['MJD-OBS']
            jd = get_jd_from_mjd(mjd)
            jds.append(jd)

        return np.array(jds)

    @property
    def nside(self):
        return list(hp.npix2nside(len(df)) for df in self.sky_probs)


    @property
    def sky_probs(self):
        """
        """
        if self._sky_probs is None:
            _x = [] 
            for (fname, threshold) in zip(self.fnames, self.prob_fractions_threshold):
                df = self.get_skyprobs(fname)
                df['in_pixels'] = 0
                idx = df.query('cumulative_prob > 1.0 - @threshold').reset_index()['hid']
                df.loc[idx, 'in_pixels'] = 1
                _x.append(df)
            self._sky_probs = _x
        return self._sky_probs

    @property
    def sky_pix(self):
        """
        """
        if self._skypix is None:
            self._skypix = list(set(df.query('in_pixels == 1').reset_index()['hid'])
                                for df in self.sky_probs)
        return self._skypix

    @staticmethod
    def check_candidate_in(ra, dec, nside, skypix):
        """
        """
        ra = np.ravel(ra)
        dec = np.ravel(dec)
        hid = hp.ang2pix(nside, ra, dec, lonlat=True) 
        val = list()
        for idx in hid:
            if idx in skypix:
                val.append('True')
            else:
                val.append('False')

        #return list(idx for idx in hid if idx in skypix)
        return val



    @staticmethod
    def get_skymap(fname):
        # Read skymap, calculate top pixels
        top_fraction = 0.90 # limit skymap top 90% region
        skymap = hp.read_map(fname)
        npix = len(skymap)
        nside = hp.npix2nside(npix)
        
        # Convert to astropy Table, easier for manipulation
        indices = np.arange(len(skymap))
        tm = Table(data=(indices, skymap), names=('hid', 'prob'))
        tm.sort('prob')
        cs = np.cumsum(tm['prob'])
        cs.name='cumulative_prob'
        tm.add_column(cs)
        
        top_pix = (tm['cumulative_prob'] > 1 - top_fraction)
        tp = Column(data=top_pix, name="top")
        tm.add_column(tp)
        top_subset = set(tm['hid'][tm['top']])
        return tm.to_pandas().set_index('hid'), top_subset
        
    @staticmethod
    def get_skyprobs(fname):
        """
        """
        data = hp.read_map(fname)
        df = pd.DataFrame(data, columns=['prob'])
        df['hid'] = np.arange(len(df))
        df.sort_values(by='prob', inplace=True, ascending=True)
        df['cumulative_prob'] = np.cumsum(df.prob)
        df.sort_values(by='hid', inplace=True)
        return df.set_index('hid')




