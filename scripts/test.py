import healpy as hp
from astropy.table import Table, Column
import numpy as np

import matplotlib.pyplot as plt
"""
Calculate if a list of objects are inside 90% localisation region of a given skymap
Varun Bhalerao, 28 April 2019

"""

# Read skymap, calculate top pixels
top_fraction = 0.90 # limit skymap top 90% region
skymap = hp.read_map("LALInference1.fits.gz,0")
npix = len(skymap)
nside = hp.npix2nside(npix)

# Convert to astropy Table, easier for manipulation
indices = np.arange(len(skymap))
tm = Table(data=(indices, skymap), names=('id', 'prob'))
tm.sort('prob')
cs = np.cumsum(tm['prob'])
cs.name='cumsum'
tm.add_column(cs)

top_pix = (tm['cumsum'] > 1 - top_fraction)
tp = Column(data=top_pix, name="top")
tm.add_column(tp)

# Cast as a set for easier comparison below
top_subset = set(tm['id'][tm['top']])
print(len(top_subset))

# Read the candidate list
hp.mollview(skymap)
plt.show()
# hp.projscatter(t['ra'], t['dec'], c='white', marker='x', lonlat=True)

#print("Candidates inside:")
#for row in t:
#    pix = hp.ang2pix(nside, row['ra'], row['dec'], lonlat=True)
#    if pix in top_subset:
#        print(row['Name'])


# print("Candidates outside:")
# for row in t:
#    pix = hp.ang2pix(nside, row['ra'], row['dec'], lonlat=True)
#    if not pix in top_subset:
#        print(row['Name'])


# Make a plot
# tnew = Table.copy(tm)
# tnew['prob'][~tnew['top']] = np.nan
# tnew.sort('id')
# hp.mollview(tnew['prob'])
# hp.projscatter(t['ra'], t['dec'], c='white', marker='x', lonlat=True)
