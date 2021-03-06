# Zach Pace, UMN/Buffalo, Summer 2013
'''USE FOR REFERENCE ONLY'''
# fits stdev, mean, and abundance amplitude for each magnitude bin
# weights bins by total bin debiased elliptical vote fraction
# fits and outputs parameters for tanhlin(s1/s2/m1/m2) and gauss(a1/a2)

def usefits(filename, ucol, rcol, votecol, (cl, ch), (rl, rh)):
	import numpy as np
	from pyfits import open
	hdulist = open(filename)
	raw = hdulist[1].data
	u, r, vote = raw.field(ucol), raw.field(rcol), raw.field(votecol)
	
	colcond = ((u - r) < ch) & ((u - r) > cl) & (r > rl) & (r < rh)
	raw = raw[colcond]
	u, r, vote = np.asfarray(raw.field(ucol), dtype = float), np.asfarray(raw.field(rcol), dtype = float), np.asfarray(raw.field(votecol), dtype = float)
	mixed_data = np.column_stack((u, r, u - r, vote))
	return mixed_data

def binme(data, mag_bins, (cl, ch), (rl, rh)):
	import numpy as np
	width = (rh - rl)/float(mag_bins)
	lowedge = np.linspace(rl, rh - width, mag_bins)
	highedge = np.linspace(rl + width, rh, mag_bins)
	edges = np.column_stack((lowedge, highedge))
	middles = 0.5 * (lowedge + highedge)

	binned_data = []
	for i in range(0, mag_bins):
		newbin = data[((data[:,1] >= edges[i,0]) & (data[:,1] < edges[i,1]))]
		binned_data.append(newbin)
	return binned_data, edges, middles

def fit_arrays(mbins):		# allows us to keep all fit values safe in an array and reference them later
	import numpy as np
	OPT = np.zeros(shape = (3, mbins, 6))	# 3d array:  6 columns, mbins rows, times 3 iterations
	ERR = np.ones(shape = (3, mbins, 6))
	return OPT, ERR

def gauss(x, s, mu, a):
	import numpy as np
	return a * np.exp(-((x - mu)**2 / (2. * s**2)))

def func0(x, s1, s2, mu1, mu2, a1, a2):
	import numpy as np
	return gauss(x, s1, mu1, a1) + gauss(x, s2, mu2, a2)

def func1(x, mu1, mu2, a1, a2):
	import numpy as np
	return func0(x, OPT[2,index,0], OPT[2,index,1], mu1, mu2, a1, a2)

def func2(x, a1, a2):
	import numpy as np
	return func1(x, OPT[2,index,2], OPT[2,index,3], a1, a2)

def g_guess(mbins):
	import numpy as np
	from sys import exit
	# g = np.array( [ [.5, .3, 2.0, 2.52, 120., 325.], [.75, .3, 2.1, 2.57, 450., 1050.], [.6, .4, 2.1, 2.55, 1800., 4800.], [.7, .4, 1.85, 2.5, 4250., 13000.], [.7, .4, 1.8, 2.5, 12000., 16500.], [.3, .2, 2.0, 2.65, 12000., 25000.], [.7, .6, 1.5, 2.5, 20000., 20000.], [.6, .8, 1.3, 2.2, 10000., 5000.], [.2, .1, 1.3, 2.3, 5000., 500.], [.8, .5, 1.1, 2.35, 5500., 2.5] ] ) # rh = -18, mbins = 10, cbins = 10
	# g = np.array( [ [.5, .3, 2.0, 2.52, 120., 325.], [.75, .3, 2.1, 2.57, 450., 1050.], [.6, .4, 2.1, 2.55, 1800., 4800.], [.7, .4, 1.85, 2.5, 4250., 13000.], [.7, .4, 1.8, 2.5, 12000., 16500.], [.3, .2, 2.0, 2.65, 12000., 25000.], [.25, .25, 1.5, 2.5, 20000., 20000.], [.3, .3, 1.5, 2.2, 12000., 5000.], [.2, .1, 1.3, 2.3, 5000., 500.] ] ) # rh = -17.7, mbins = 9, cbins = 20
	g = np.array( [ [.2, .15, 1.8, 2.55, 6., 147.], [.6, .15, 2.1, 2.55, 65., 725.], [.6, .15, 2.1, 2.55, 260., 2900.], [.6, .15, 2.1, 2.55, 1000., 8000.], [.6, .3, 1.6, 2.55, 1500., 9500.], [.3, .3, 1.4, 2.5, 2000., 7000.], [.6, .3, 1.5, 2.45, 1250., 1750.], [.2, .5, 1.3, 2., 650., 200.] ] )
	if g.shape == (mbins, 6):
		return g
	else:
		print 'Wrong dimensions!  Guess array must have', mbins, 'rows and 6 columns!'
		exit()

def t_guess():
	import numpy as np
		# g = np.array( [ [ [ 0.298, 0.014, -0.067, -19.90, 0.58 ], [ 0.152, 0.008, 0.044, -19.91, 0.94 ] ], [ [ 2.279, -0.037, -0.108, -19.81, 0.96 ], [ 1.790, -0.053, -0.363, -20.75, 1.12 ] ] ] ) # baldry's values
	g = np.array( [ [ [ 0.298, 0.014, -0.067, -19.90, 0.58 ], [ 0.152, 0.008, 0.044, -19.91, 0.94 ] ], [ [ 1.790, -0.053, -0.363, -20.75, 1.12 ] , [ 2.279, -0.037, -0.108, -19.81, 0.96 ] ] ] )
	return g

def fitmagbin(data, mbin, cbins, iteration):
	import numpy as np
	from scipy.optimize import curve_fit

	N, B = np.histogram(data[ : , 2], bins = cbins, weights = data[ : , 3])
	# print np.sum(np.isnan(N)) # no nans show up
	E = np.sqrt(N)
	E = np.where(E > 0., E, 2.)
	# print E
	bincenters = 0.5 * (B[1:] + B[:-1])

	if iteration == 0:
		opt, cov = curve_fit(func0, bincenters, N, p0 = g_guess[mbin, 0: ], sigma = E, maxfev = 100000)
	if iteration == 1:
		opt, cov = curve_fit(func1, bincenters, N, p0 = g_guess[mbin, 2: ], sigma = E, maxfev = 100000)
	if iteration == 2:
		opt, cov = curve_fit(func2, bincenters, N, p0 = g_guess[mbin, 4: ], sigma = E, maxfev = 100000)

	err = np.sqrt(np.diagonal(cov))

	ERR[iteration, mbin, 2 * iteration : ] = err
	OPT[iteration, mbin, 2 * iteration : ] = opt

	print 'fit values'
	fitdata = np.column_stack((opt, err))
	print fitdata
	return N, E, bincenters

def plotmagbin(data, best_fit, iteration, mbin, E, bincenters):
	import numpy as np
	import matplotlib.pyplot as plt

	x = np.linspace(cl, ch, num = 200)
	gauss1 = gauss(x, best_fit[iteration, mbin, 0], best_fit[iteration, mbin, 2], best_fit[iteration, mbin, 4])
	gauss2 = gauss(x, best_fit[iteration, mbin, 1], best_fit[iteration, mbin, 3], best_fit[iteration, mbin, 5])
	plt.plot(x, gauss1, 'b--', x, gauss2, 'r--', x, gauss1 + gauss2, 'g')
	
	'''
	ggauss1 = gauss(x, g_guess[mbin, 0], g_guess[mbin, 2], g_guess[mbin, 4])
	ggauss2 = gauss(x, g_guess[mbin, 1], g_guess[mbin, 3], g_guess[mbin, 5])
	plt.plot(x, ggauss1, 'b:', x, ggauss2, 'r:', x, ggauss1 + ggauss2, 'g:')
	'''

	plt.errorbar(bincenters, data, yerr = E, marker = 'o', fmt = 'x', color = 'g', ecolor = 'g', capsize = 6)
	plt.xlabel('color (u - r)')
	plt.ylabel('frequency')
	plt.title('iteration %s, bin %s (%s < r < %s)' %(iteration, index, edges[index, 0], edges[index, 1]))
	plt.grid(True)
	# plt.legend(('blue cloud', 'red series', 'envelope gaussian'), loc = 'best')
	plt.show()

def tanhlin(x, p0, p1, q0, q1, q2):
	import numpy as np
	return p0 + p1 * (x + 20.) + q0 * np.tanh((x - q1)/q2)

def tanhlinfit(blue_col, red_col):
	from scipy.optimize import curve_fit
	bluefit, bluefit_e = curve_fit(tanhlin, middles, OPT[iteration, : , blue_col], sigma = ERR[iteration, : , blue_col], p0 = t_guess[iteration, 0], maxfev = 100000)
	redfit, redfit_e = curve_fit(tanhlin, middles, OPT[iteration, : ,  red_col], sigma = ERR[iteration, : , red_col], p0 = t_guess[iteration, 0], maxfev = 100000)
	return bluefit, redfit

def tanhlinplot(bluefit, redfit, blue_col, red_col):
	import matplotlib.pyplot as plt

	m = np.linspace(rl, rh, 200)
	plt.plot(m, tanhlin(m, bluefit[0], bluefit[1], bluefit[2], bluefit[3], bluefit[4]), color = 'b', linestyle = '-')
	plt.plot(m, tanhlin(m, redfit[0], redfit[1], redfit[2], redfit[3], redfit[4]), color = 'r', linestyle = '-')

	gtan1 = tanhlin(m, t_guess[iteration, 0, 0], t_guess[iteration, 0, 1], t_guess[iteration, 0, 2], t_guess[iteration, 0, 3], t_guess[iteration, 0, 4])
	gtan2 = tanhlin(m, t_guess[iteration, 1, 0], t_guess[iteration, 1, 1], t_guess[iteration, 1, 2], t_guess[iteration, 1, 3], t_guess[iteration, 1, 4])
	plt.plot(m, gtan1, 'b:', m, gtan2, 'r:')

	#	first plot the blue cloud data points
	plt.errorbar(middles, OPT[iteration, : , blue_col], yerr = ERR[iteration, : , blue_col], marker = 'd', c = 'b', fmt = 'd', ecolor = 'b')
	#	now plot the red cloud data points
	plt.errorbar(middles, OPT[iteration, : , red_col], yerr = ERR[iteration, : , red_col], marker = 'x', c = 'r', fmt = 'x', ecolor = 'r')
	

def plotoptions(xlabel, ylabel, title, grid, legend, legendtext):
	import matplotlib.pyplot as plt
	
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.grid(grid)
	if legend == True:
		plt.legend(legendtext, loc = 'best')
	plt.show()

#	=====
#	call functions
#	=====

import numpy as np

cl = 0.
ch = 3.5
rl = -24.
rh = -17.7
mbins = 8
cbins = 20
votecol = 4

mixed_data = usefits('../../Downloads/Zoo1_CM2_zjpace.fit', 6, 7, 4, (cl, ch), (rl, rh))	# returns array( [ u, r, u - r] )

binned_data, edges, middles = binme(mixed_data, mbins, (cl, ch), (rl, rh))	# returns list of arrays, binned by r magnitude
print edges, middles

OPT, ERR = fit_arrays(mbins)	# returns an empty 3-deep, mbins-down, 6-across array for each the optimized parameters and the error matrix

g_guess = g_guess(mbins)

# print g_guess

# =====
# iteration 0 (fit all, fix sigma)
# =====

iteration = 0
for index, b in enumerate(binned_data):
	print 'color bin', index
	N, E, bincenters = fitmagbin(b, index, cbins, iteration)
	# print OPT[iteration, : , 0]	# watch the OPT array being built
	plotmagbin(N, OPT, iteration, index, E, bincenters)

t_guess = t_guess()

bluefit_sig, redfit_sig = tanhlinfit(0, 1)
'''print 'blue:'
print bluefit_sig
print 'red:'
print redfit_sig'''

S = np.array( [ tanhlin(middles, bluefit_sig[0], bluefit_sig[1],  bluefit_sig[2],  bluefit_sig[3],  bluefit_sig[4]) , tanhlin(middles, redfit_sig[0], redfit_sig[1],  redfit_sig[2],  redfit_sig[3],  redfit_sig[4]) ] )

# print S
# print OPT[0, : , 0]

OPT[2, : , 0] = S[0, : ]
OPT[2, : , 1] = S[1, : ]
OPT[1, : , 0] = S[0, : ]
OPT[1, : , 1] = S[1, : ]

ERR[2, : , 0] = ERR[0, : , 0]
ERR[2, : , 1] = ERR[0, : , 1]
ERR[1, : , 0] = ERR[0, : , 0]
ERR[1, : , 1] = ERR[0, : , 1]

# print 'OPT:', '\n', OPT
# print 'ERR:', '\n', ERR

tanhlinplot(bluefit_sig, redfit_sig, 0, 1)
plotoptions('r magnitude', 'st. dev.', 'color spread of magnitude bins', True, True, ('blue series computer fit', 'red series computer fit', 'blue series Baldry fit', 'red series Baldry fit'))

# =====
# iteration 1 (fit mean & amplitude, fix mean)
# =====

iteration = 1
for index, b in enumerate(binned_data):
	print 'color bin', index
	N, E, bincenters = fitmagbin(b, index, cbins, iteration)
	# plotmagbin(N, OPT, iteration, index, E, bincenters)
	
bluefit_mu, redfit_mu = tanhlinfit(2,3)

'''print 'blue:'
print bluefit_mu
print 'red:'
print redfit_mu'''

M = np.array( [ tanhlin(middles, bluefit_mu[0], bluefit_mu[1],  bluefit_mu[2],  bluefit_mu[3],  bluefit_mu[4]) , tanhlin(middles, redfit_mu[0], redfit_mu[1],  redfit_mu[2],  redfit_mu[3],  redfit_mu[4]) ] )

OPT[2, : , 2] = M[0, : ]
OPT[2, : , 3] = M[1, : ]

ERR[2, : , 2] = ERR[1, : , 2]
ERR[2, : , 3] = ERR[1, : , 3]

tanhlinplot(bluefit_mu, redfit_mu, 2, 3)
plotoptions('r magnitude', 'mean color (u - r)', 'mean color of magnitude bins', True, True, ('blue series computer fit', 'red series computer fit', 'blue series Baldry fit', 'red series Baldry fit'))

# =====
# iteration 2 (fit amplitude, fix amplitude)
# =====

iteration = 2
for index, b in enumerate(binned_data):
	print 'color bin', index
	N, E, bincenters = fitmagbin(b, index, cbins, iteration)
	# plotmagbin(N, OPT, iteration, index, E, bincenters)

	#	first plot the blue cloud data points
import matplotlib.pyplot as plt

final_params = OPT[2]
final_errors = ERR[2]

import numpy as np
from scipy.optimize import curve_fit
bluefit_a, bluefit_a_e = curve_fit(gauss, middles, OPT[2, : , 4], p0 = (1., -19.6, 13500.) )
redfit_a, redfit_a_e = curve_fit(gauss, middles, OPT[2, : , 5], p0 = (1., -20.5, 18500.) )

'''print 'blue:'
print bluefit_a
print 'red:'
print redfit_a'''

x = np.linspace(rl, rh, 200)
plt.plot(x, gauss(x, bluefit_a[0], bluefit_a[1], bluefit_a[2]), c = 'b')
plt.plot(x, gauss(x, redfit_a[0], redfit_a[1], redfit_a[2]), c = 'r')

plt.errorbar(middles, OPT[iteration, : , 4], yerr = ERR[iteration, : , 4], marker = 'd', c = 'b', fmt = 'd', ecolor = 'b')
	#	now plot the red cloud data points
plt.errorbar(middles, OPT[iteration, : , 5], yerr = ERR[iteration, : , 5], marker = 'x', c = 'r', fmt = 'x', ecolor = 'r')

plotoptions('r magnitude', 'counts at peak', 'color abundance', True, True, ('blue fit', 'red fit'))

A = np.array( [ gauss(middles, bluefit_a[0], bluefit_a[1], bluefit_a[2]), gauss(middles, redfit_a[0], redfit_a[1], redfit_a[2]) ] )

print A

final_params[ : , 4] = A[0]
final_params[ : , 5] = A[1]

# =====
# Now try plotting a 3d double-gaussian
# =====

mag1d = np.linspace(rl, rh, mbins)
c1d = np.linspace(cl, rh, cbins)
print '\n', 'FINAL PARAMETER FITS:', '\n'
print bluefit_sig
print bluefit_mu
print bluefit_a
print redfit_sig
print redfit_mu
print redfit_a

r = np.linspace(rl, rh, 200)
bluemean = tanhlin(r, bluefit_mu[0], bluefit_mu[1], bluefit_mu[2], bluefit_mu[3], bluefit_mu[4])
redmean = tanhlin(r, redfit_mu[0], redfit_mu[1], redfit_mu[2], redfit_mu[3], redfit_mu[4])
bluesig = tanhlin(r, bluefit_sig[0], bluefit_sig[1], bluefit_sig[2], bluefit_sig[3], bluefit_sig[4])
redsig = tanhlin(r, redfit_sig[0], redfit_sig[1], redfit_sig[2], redfit_sig[3], redfit_sig[4])
plt.xlabel('r magnitude')
plt.ylabel('color (u - r)')
plt.title('means of color distributions, +/- two-sigma')
plt.plot(r, bluemean - 2.*bluesig, 'b:', r, bluemean - bluesig, 'b--', r, bluemean, 'b-', r, bluemean + bluesig, 'b--', r, bluemean + 2.*bluesig, 'b:', r, redmean - 2.*redsig, 'r:', r, redmean - redsig, 'r--', r, redmean, 'r-', r, redmean + redsig, 'r--', r, redmean + 2.*redsig, 'r:')
plt.show()

def f_r(r, c):
	import numpy as np
	a = gauss(r, redfit_a[0], redfit_a[1], redfit_a[2])
	m = tanhlin(r, redfit_mu[0], redfit_mu[1], redfit_mu[2], redfit_mu[3], redfit_mu[4])
	s = tanhlin(r, redfit_sig[0], redfit_sig[1], redfit_sig[2], redfit_sig[3], redfit_sig[4])
	return gauss(c, s, m, a)

def f_b(r, c):
	import numpy as np
	a = gauss(r, bluefit_a[0], bluefit_a[1], bluefit_a[2])
	m = tanhlin(r, bluefit_mu[0], bluefit_mu[1], bluefit_mu[2], bluefit_mu[3], bluefit_mu[4])
	s = tanhlin(r, bluefit_sig[0], bluefit_sig[1], bluefit_sig[2], bluefit_sig[3], bluefit_sig[4])
	return gauss(c, s, m, a)

RBs = 250
CBs = 250

r_raw = np.linspace(rl, -16., RBs)
c_raw = np.linspace(cl, ch, CBs)

r_array = np.tile(r_raw, (CBs, 1))
c_array = np.tile(c_raw, (RBs, 1)).transpose()

Z_r = f_r(r_array, c_array)
Z_b = f_b(r_array, c_array)

Z_tot = Z_r + Z_b

from mpl_toolkits.mplot3d import axes3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(r_array, c_array, Z_tot/np.amax(Z_tot), rstride = 8, cstride = 8, alpha = .75, cmap = 'spectral')
plt.show()
