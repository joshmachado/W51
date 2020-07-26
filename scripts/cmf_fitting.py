from astropy.table import Table
from astropy import units as units
import powerlaw
import pylab as pl
import numpy as np
from astropy import stats
import astropy


fp = '/Users/josh/GitHub/W51/'
t = Table.read(fp+'data/dendro_catalog.tex')
masses = t['mass']
masses = masses[1:len(masses)]
masses = np.asarray(masses).astype('float64')
masses = np.trim_zeros(np.sort(masses))
masses = masses[~np.isnan(masses)]

error = t['mass_uncertainty']
#error = error[1:len(error)]
error = np.asarray(error).astype('float64')
error = error[1:len(error)]
error = np.trim_zeros(np.sort(error))
error = error[~np.isnan(error)]


fit = powerlaw.Fit(masses)
alpha = fit.power_law.alpha


#determining optimal bin size
blocks = astropy.stats.histogram(masses, bins='blocks')

pl.hist(masses, log=True, bins=np.logspace(np.log10(0.1),np.log10(np.amax(masses)), 8))
pl.xscale('log')
pl.yscale('log')
pl.show()


#CMF
pl.plot(np.log10(np.sort(masses)),np.log10(np.linspace(1,0,len(masses), endpoint=False)))
#Fit lines to cmf
pl.show()



#Estimate error on alpha
percent_error = error / masses
num_trials = int(1e4)
alpha_est = []
for i in range(num_trials):
	completemass_werror = masses + (percent_error*masses) * np.random.randn(len(masses))
	#completemass_werror = completemass_werror[completemass_werror > 1.6]
	fit = powerlaw.Fit(completemass_werror)
	alpha_est.append(fit.alpha)

bootstrapped_alpha = np.mean(alpha_est)
bootstrapped_sigma = np.std(alpha_est)