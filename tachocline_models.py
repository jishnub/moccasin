import numpy as np
# from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import units
import itertools

sample_size = 10

center_RH = np.array([0.692,0.705,0.6947,0.695,0.691,0.693,0.697,0.6916])
errors_center_RH = np.array([0.005,0.0027,0.0035,0.005,0.004,0.002,0.002,0.0019])

centers_sampled = np.zeros([center_RH.size*sample_size])

width_RH = np.array([0.09,0.048,0.033,0.05,0.01,0.039,0.019,0.0162])
errors_width_RH = np.array([0.04,0.0127,0.0069,0.03,0.03,0.002,0.001,0.0032])

width_sampled = np.zeros([width_RH.size*sample_size])

for result_no,(c,err) in enumerate(zip(center_RH,errors_center_RH)):
	centers_sampled[result_no*sample_size:(result_no+1)*sample_size] = np.random.normal(c,err,sample_size)

for result_no,(w,err) in enumerate(zip(width_RH,errors_width_RH)):
	width_sampled[result_no*sample_size:(result_no+1)*sample_size] = np.random.normal(w,err,sample_size)

cov_center_width = np.cov(centers_sampled,width_sampled)

realizations = np.random.multivariate_normal([np.mean(centers_sampled),np.mean(width_sampled)],
											cov_center_width,1000)

centers,widths = realizations.T
nmodels = 100
realizations = realizations[(0.675<centers) & (centers<0.725) & (0.005<widths) & (widths<0.05)][:nmodels]
centers,widths = realizations.T

nmodels = centers.size


radius = fits.open("radius.fits")[0].data.squeeze()
nr = radius.size

omega_inner = 440*units.nHz
omega_outer_equator = 454*units.nHz
domega = omega_outer_equator - omega_inner

smin = 1; smax = 1; ns = (smax-smin)//2 + 1

wmodels = np.zeros((nmodels,ns,nr))*units.uHz

for modelno,(c,w) in enumerate(realizations[:nmodels]):
	# W1 is radial differential rotation
	wmodels[modelno,0] = omega_inner + domega*(np.tanh((radius - c)/w)+1)/2

fits.writeto("wmodels.fits",wmodels.value,overwrite=True)