from astropy.table import Table
from astropy.io import ascii, fits
from astropy import wcs
from astropy import coordinates
import matplotlib.pylab as plt
cat_w51e = ascii.read('/Users/josh/GitHub/W51/data/W51-E.sw.sources.fin.ok_dec2020.cat', data_start=0, format='commented_header', header_start=110,  comment="!")
cat_w51irs2 = ascii.read('/Users/josh/GitHub/W51/data/W51-IRS2.sw.sources.fin.ok_dec2020.cat', data_start=0, format='commented_header', header_start=110, comment="!")
parmap = fits.open('/Users/josh/GitHub/W51/data/par_maps.fits')
ww = wcs.WCS(parmap[0].header).celestial
plt.figure(figsize=(20,20))
plt.subplot(projection=ww)
ax = plt.gca()
plt.imshow(parmap[0].data[0])
plt.plot(cat_w51e['WCS_ACOOR'], cat_w51e['WCS_DCOOR'], 'w.', transform=ax.get_transform('world'))
plt.plot(cat_w51irs2['WCS_ACOOR'], cat_w51irs2['WCS_DCOOR'], 'r.', transform=ax.get_transform('world'))
xpix, ypix = ww.wcs_world2pix(cat_w51e['WCS_ACOOR'], cat_w51e['WCS_DCOOR'], 0)
temperatures_w51e = parmap[0].data[0][ypix.astype('int'), xpix.astype('int')]
xpix, ypix = ww.wcs_world2pix(cat_w51irs2['WCS_ACOOR'], cat_w51irs2['WCS_DCOOR'], 0)
temperatures_irs2 = parmap[0].data[0][ypix.astype('int'), xpix.astype('int')]
