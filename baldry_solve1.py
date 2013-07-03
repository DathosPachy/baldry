# Zach Pace, UMN/Buffalo, Summer 2013

# fits stdev, mean, and abundance amplitude for each magnitude bin
# weights bins by total bin based on either flag value OR debiased vote fraction, based on votecol specified in usefits
# fits and outputs parameters for tanhlin(s1/m1)

def usefits(filename, ucol, rcol, votecol, (cl, ch), (rl, rh)):
	import numpy as np
	from pyfits import open
	hdulist = open(filename)
	raw = hdulist[1].data
	u, r = raw.field(ucol), raw.field(rcol)
	cols = hdulist[1].columns
	colsnames = cols.names
	print 'Selecting based on column', votecol, ':', colsnames[votecol]
	
	colcond = ((u - r) < ch) & ((u - r) > cl) & (r > rl) & (r < rh)
	raw = raw[colcond]
	u, r, vote = np.asfarray(raw.field(ucol), dtype = float), np.asfarray(raw.field(rcol), dtype = float), np.asfarray(raw.field(votecol), dtype = float)
	mixed_data = np.column_stack((u, r, u - r, vote)) # use for vote morphology
	# mixed_data = np.column_stack((u, r, u - r)) # use for base case
	# print 'mixed data:', '\n', mixed_data
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
	OPT = np.zeros(shape = (3, mbins, 3))	# 3d array:  6 columns, mbins rows, times 3 iterations
	ERR = np.ones(shape = (3, mbins, 3))
	return OPT, ERR

def gauss(x, s, mu, a):
	import numpy as np
	return a * np.exp(-((x - mu)**2 / (2. * s**2)))

def gauss1(x, mu, a):
	import numpy as np
	return gauss(x, OPT[2, index, 0], mu, a)

def gauss2(x, a):
	import numpy as np
	return gauss1(x, OPT[2, index, 1], a)

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
	g = np.array( [ [.2, 2.5, 140.], [.15, 2.55, 650.], [.15, 2.55, 2350.], [.15, 2.5, 3000.], [.6, 2.5, 1500.], [.3, 2.4, 525.], [.6, 2.35, 120.], [.2, 2.4, 14.] ] ) # USE FOR ELLIPTICAL FLAG (votecol = 1)
	if g.shape == (mbins, 3):
		return g
	else:
		print 'Wrong dimensions!  Guess array must have', mbins, 'rows and 3 columns!'
		exit()

def t_guess():
	import numpy as np
	# g = np.array( [ [ [ 0.298, 0.014, -0.067, -19.90, 0.58 ], [ 0.152, 0.008, 0.044, -19.91, 0.94 ] ], [ [ 2.279, -0.037, -0.108, -19.81, 0.96 ], [ 1.790, -0.053, -0.363, -20.75, 1.12 ] ] ] ) # baldry's values
	g = np.array( [ [ [ 0.298, 0.014, -0.067, -19.90, 0.58 ], [ 0.152, 0.008, 0.044, -19.91, 0.94 ] ], [ [ 1.790, -0.053, -0.363, -20.75, 1.12 ] , [ 2.279, -0.037, -0.108, -19.81, 0.96 ] ] ] )
	return g

def fitmagbin(data, mbin, cbins, iteration):
	import numpy as np
	from scipy.optimize import curve_fit
	# N, B = np.histogram(data[ : , 2], bins = cbins) # use for base case 
	N, B = np.histogram(data[ : , 2], bins = cbins, weights = data[ : , 3]) # use for morphology
	# print np.sum(np.isnan(N)) # no nans show up
	E = np.sqrt(N)
	E = np.where(E > 0., E, 2.)
	# print E
	bincenters = 0.5 * (B[1:] + B[:-1])

	if iteration == 0:
		opt, cov = curve_fit(gauss, bincenters, N, p0 = g_guess[mbin, 0: ], sigma = E, maxfev = 100000)
	if iteration == 1:
		opt, cov = curve_fit(gauss1, bincenters, N, p0 = g_guess[mbin, 1: ], sigma = E, maxfev = 100000)
	if iteration == 2:
		opt, cov = curve_fit(gauss2, bincenters, N, p0 = g_guess[mbin, 2: ], sigma = E, maxfev = 100000)

	err = np.sqrt(np.diagonal(cov))
	err = np.absolute(err)
	ERR[iteration, mbin, iteration : ] = err
	opt = np.absolute(opt)
	OPT[iteration, mbin, iteration : ] = opt

	print 'fit values'
	fitdata = np.column_stack((opt, err))
	print fitdata
	return N, E, bincenters

def plotoptions(xlabel, ylabel, title, grid, legend, legendtext):
	import matplotlib.pyplot as plt
	
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.grid(grid)
	if legend == True:
		plt.legend(legendtext, loc = 'best')
	# plt.savefig(title + '.png')
	plt.show()

def plotmagbin(data, best_fit, iteration, mbin, E, bincenters):
	import numpy as np
	import matplotlib.pyplot as plt

	x = np.linspace(cl, ch, num = 200)
	gauss1 = gauss(x, best_fit[iteration, mbin, 0], best_fit[iteration, mbin, 1], best_fit[iteration, mbin, 2])

	plt.plot(x, gauss1, 'k--')
	
	'''
	ggauss1 = gauss(x, g_guess[mbin, 0], g_guess[mbin, 2], g_guess[mbin, 4])
	ggauss2 = gauss(x, g_guess[mbin, 1], g_guess[mbin, 3], g_guess[mbin, 5])
	plt.plot(x, ggauss1, 'b:', x, ggauss2, 'r:', x, ggauss1 + ggauss2, 'g:')
	'''

	plt.errorbar(bincenters, data, yerr = E, marker = 'o', fmt = 'x', color = 'g', ecolor = 'g', capsize = 6)
	plotoptions('color (u - r)', 'frequency', 'iteration %s, bin %s (%s < r < %s)' %(iteration, index, edges[index, 0], edges[index, 1]), True, False, ('envelope gaussian'))

def tanhlin(x, p0, p1, q0, q1, q2):
	import numpy as np
	return p0 + p1 * (x + 20.) + q0 * np.tanh((x - q1)/q2)

def tanhlinfit(blue_col):
	from scipy.optimize import curve_fit
	bluefit, bluefit_e = curve_fit(tanhlin, middles[:-1], OPT[iteration, :-1 , blue_col], sigma = ERR[iteration, : -1, blue_col], p0 = t_guess[iteration, 0], maxfev = 100000)

	return bluefit

def tanhlinplot(bluefit, blue_col):
	import matplotlib.pyplot as plt

	m = np.linspace(rl, rh, 200)
	plt.plot(m, tanhlin(m, bluefit[0], bluefit[1], bluefit[2], bluefit[3], bluefit[4]), color = 'b', linestyle = '-')

	#	first plot the blue cloud data points
	plt.errorbar(middles, OPT[iteration, : , blue_col], yerr = ERR[iteration, : , blue_col], marker = 'd', c = 'b', fmt = 'd', ecolor = 'b')

def makelogfile(zooversion, destination, target, method):	# only writes to a subdirectory of where the script is originally running!
	# target: 'what population are we supposed to be targeting?' (STRING: ell, spi, or all)
	# method: 'are we using flags or debiased count?' (STRING: flag, deb, or none)
	import time, os, pickle, numpy as np
	timestr = time.strftime('%Y-%m-%d-%H-%M')
	filename = 'z' + str(zooversion) + '-' + target + '-' + method + '-' + timestr + '.txt'
	filename = os.path.join(destination, filename)
	print 'Creating log file: ', filename
	f = open(filename, 'wb')
	
	f.write(timestr + '\n')
	f.write(str(rl) + '\n')
	f.write(str(rh) + '\n')
	f.write(str(cl) + '\n')
	f.write(str(ch) + '\n')
	f.write(str(mbins) + '\n')
	f.write(str(cbins) + '\n')
	'''f.write(str(votecol) + '\n')'''
	f.write(str(bluefit_sig) + '\n')
	f.write(str(redfit_sig) + '\n')
	f.write(str(bluefit_mu) + '\n')
	f.write(str(redfit_mu) + '\n')
	f.write(str(bluefit_a) + '\n')
	f.write(str(redfit_a) + '\n')

	comment = str(raw_input('Enter any other comments now: '))
	f.write(comment)

	print 'Log file done! Accessible at', filename

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
votecol = 1

mixed_data = usefits('../../Downloads/Zoo1_CM2_zjpace.fit', 6, 7, votecol, (cl, ch), (rl, rh))	# returns array( [ u, r, u - r] )

samplesize = len(mixed_data)
print 'Size of Sample:', samplesize

binned_data, edges, middles = binme(mixed_data, mbins, (cl, ch), (rl, rh))	# returns list of arrays, binned by r magnitude
print edges, middles

OPT, ERR = fit_arrays(mbins)	# returns an empty 3-deep, mbins-down, 3-across array for each the optimized parameters and the error matrix

g_guess = g_guess(mbins)

# print g_guess

# =====
# iteration 0 (fit all, fix sigma)
# =====

iteration = 0
for index, b in enumerate(binned_data):
	print 'color bin', index
	N, E, bincenters = fitmagbin(b, index, cbins, iteration)
	# plotmagbin(N, OPT, iteration, index, E, bincenters)

t_guess = t_guess()

bluefit_sig = tanhlinfit(0)
'''print 'blue:'
print bluefit_sig
print 'red:'
print redfit_sig'''

S = np.array( [ tanhlin(middles, bluefit_sig[0], bluefit_sig[1],  bluefit_sig[2],  bluefit_sig[3],  bluefit_sig[4]) ] )

# print S
# print OPT[0, : , 0]

OPT[2, : , 0] = S
OPT[1, : , 0] = S

ERR[2, : , 0] = ERR[0, : , 0]
ERR[1, : , 0] = ERR[0, : , 0]

# print 'OPT:', '\n', OPT
# print 'ERR:', '\n', ERR

tanhlinplot(bluefit_sig, 0)
plotoptions('r magnitude', 'st. dev.', 'color spread of magnitude bins', True, True, ('combined series fit'))

# =====
# iteration 1 (fit mean & amplitude, fix mean)
# =====

iteration = 1
for index, b in enumerate(binned_data):
	print 'color bin', index
	N, E, bincenters = fitmagbin(b, index, cbins, iteration)
	# plotmagbin(N, OPT, iteration, index, E, bincenters)
	
bluefit_mu = tanhlinfit(1)

'''print 'blue:'
print bluefit_mu
print 'red:'
print redfit_mu'''

M = np.array( [ tanhlin(middles, bluefit_mu[0], bluefit_mu[1],  bluefit_mu[2],  bluefit_mu[3],  bluefit_mu[4]) ] )

OPT[2, : , 2] = M[0, : ]

ERR[2, : , 2] = ERR[1, : , 2]

tanhlinplot(bluefit_mu, 1)
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
bluefit_a = curve_fit(gauss, middles, OPT[2, : , 2], p0 = (1., -19.6, 13500.) )

'''print 'blue:'
print bluefit_a
print 'red:'
print redfit_a'''

print bluefit_a

x = np.linspace(rl, rh, 200)
plt.plot(x, gauss(x, bluefit_a[0], bluefit_a[1], bluefit_a[2]), c = 'b')

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

'''def clims(r):
	return tanhlin(r, bluefit_mu[0], bluefit_mu[1], bluefit_mu[2], bluefit_mu[3], bluefit_mu[4]), tanhlin(r, redfit_mu[0], redfit_mu[1], redfit_mu[2], redfit_mu[3], redfit_mu[4])

def min_at_mag(r):
	import numpy as np
	color = np.linspace(clims(r)[0], clims(r)[1], 100)
	vals = f_r(r, color) + f_b(r, color)
	mag_min = np.min(vals)
	i = list(vals).index(mag_min)
	min_color = color[i]
	return min_color

def fit_mtn_pass(r, magsamples):
	import numpy as np
	from scipy.optimize import curve_fit
	samplemags = np.linspace(rl, rh, magsamples)
	mins_at_mags = min_at_mag(samplemags)
	opt = curve_fit(tanhlin, r, mins_at_mags, p0 = 0.5 * (t_guess[0] + t_guess[1]), maxfev = 100000)
	return opt

samples = 10

r_test = np.linspace(-24., -17.7, samples)

split = fit_mtn_pass(r_test, samples)'''

plt.figure()
m = np.linspace(rl, rh, 200)
plt.plot(m, tanhlin(m, bluefit_mu[0], bluefit_mu[1], bluefit_mu[2], bluefit_mu[3], bluefit_mu[4]), color = 'b', linestyle = '-')
plt.plot(m, tanhlin(m, redfit_mu[0], redfit_mu[1], redfit_mu[2], redfit_mu[3], redfit_mu[4]), color = 'r', linestyle = '-')
plt.plot(m, 0.5 * (tanhlin(m, bluefit_mu[0], bluefit_mu[1], bluefit_mu[2], bluefit_mu[3], bluefit_mu[4]) + tanhlin(m, redfit_mu[0], redfit_mu[1], redfit_mu[2], redfit_mu[3], redfit_mu[4])), color = 'g', linestyle = '-')

# im = plt.imshow(Z_r + Z_b, origin = 'lower')
im = plt.pcolor(r_array, c_array, Z_r + Z_b, cmap = 'terrain')
CS = plt.contour(r_array, c_array, Z_r + Z_b, 10, colors = 'k')
plt.clabel(CS, inline = 1, fmt = '%1.0f', fontsize = 8, colors = 'k')	# only for contour/contourf
CB = plt.colorbar(CS, shrink = 0.8, extend = 'both')
CBI = plt.colorbar(im, orientation = 'horizontal', shrink = 0.8)

l,b,w,h = plt.gca().get_position().bounds
ll,bb,ww,hh = CB.ax.get_position().bounds
CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])
plotoptions('r magnitude', 'color', 'Distribution of galaxies', True, True, ('mean color of red series','mean color of blue cloud','split function'))
plt.show()

makelogfile(1, 'zoo1/flags/spi', 'spi', 'flags')
