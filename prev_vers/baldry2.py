# Zach Pace, UMN/Buffalo, Summer 2013
# takes elliptical data and fits some curves to the CMD, based on raw data

# NOTE: uses 10 bins; to do anything different requires modification of code

def tanhlin(x, p0, p1, q0, q1, q2):
	import numpy as np
	return p0 + p1 * (x + 20.) + q0 * np.tanh((x - q1)/q2)

def gauss(x, s, mu, a):
	import numpy as np
	return a * np.exp(-((x - mu)**2 / (2. * s**2)))

def func0(x, s1, s2, mu1, mu2, a1, a2):
	import numpy as np
	return a1 * np.exp(-((x - mu1)**2 / (2. * s1**2))) + a2 * np.exp(-((x - mu2)**2 / (2. * s2**2)))

import numpy as np	# this bit is actually necessary, it appears
sigma1_f = np.array( [ ] )
sigma2_f = np.array( [ ] )

def func1_n(x, n, mu1, mu2, a1, a2):
	import numpy as np
	return a1 * np.exp(-((x - mu1)**2 / (2. * sigma1_f[n]**2))) + a2 * np.exp(-((x - mu2)**2 / (2. * sigma2_f[n]**2)))

def cmfit(filename, ucol, rcol, c_low, c_high, r_low, r_high):
	# dependencies
	import pyfits as pf
	import numpy as np
	from scipy.optimize import curve_fit
	from scipy.stats import chisquare
	import matplotlib.pyplot as plt
	# getting file, and extracting data
	hdulist = pf.open(filename)
	data = hdulist[1].data
	u = data.field(ucol)
	r = data.field(rcol)

	colcond = ((u - r) < c_high) & ((u - r) > c_low) & (r > r_low) & (r < r_high)
	data = data[colcond]
	u = data.field(ucol)
	u = np.asfarray(u,dtype=float)
	r = data.field(rcol)
	r = np.asfarray(r,dtype=float)
	color = u - r
	cmdata = np.column_stack((u, r, color))

	# -----
	# bin data (since none of the histogram functions preserve non-binned data correctly)
	# -----

	edges = np.linspace(r_low, r_high, 10+1) # 10 bins
	bc0 = (r >= edges[0]) & (r < edges[1])
	bc1 = (r >= edges[1]) & (r < edges[2])
	bc2 = (r >= edges[2]) & (r < edges[3])
	bc3 = (r >= edges[3]) & (r < edges[4])
	bc4 = (r >= edges[4]) & (r < edges[5])
	bc5 = (r >= edges[5]) & (r < edges[6])
	bc6 = (r >= edges[6]) & (r < edges[7])
	bc7 = (r >= edges[7]) & (r < edges[8])
	bc8 = (r >= edges[8]) & (r < edges[9])
	bc9 = (r >= edges[9]) & (r <= edges[10])

	b0_0 = cmdata[bc0]	# bin zero, iteration 0
	b1_0 = cmdata[bc1]
	b2_0 = cmdata[bc2]
	b3_0 = cmdata[bc3]
	b4_0 = cmdata[bc4]
	b5_0 = cmdata[bc5]
	b6_0 = cmdata[bc6]
	b7_0 = cmdata[bc7]
	b8_0 = cmdata[bc8]
	b9_0 = cmdata[bc9]

	# -----
	# now create a histogram for each bin
	# -----

	b_all = [ b0_0, b1_0, b2_0, b3_0, b4_0, b5_0, b6_0, b7_0, b8_0, b9_0 ]
	p0 = np.array( [ [.2, .15, 1.8, 2.55, 80., 245.], [.6, .15, 2.1, 2.55, 325., 925.], [.6, .15, 2.1, 2.55, 325., 925.], [.6, .15, 2.1, 2.55, 325., 925.], [.6, .3, 1.6, 2.55, 15000., 22000.], [.3, .3, 1.4, 2.5, 10000., 8000.], [.6, .3, 1.5, 2.45, 9000., 10000.], [.2, .5, 1.3, 2., 1700., 90.], [.2, .1, 1.15, 2.1, 1300., 10.], [.1, .3, 1.1, 2.35, 60., 2.5] ] )

	sigma1_0 = np.array( [ ] )
	sigma2_0 = np.array( [ ] )
	mag = np.array( [ ] )
	sigma1e_0 = np.array( [ ] )
	sigma2e_0 = np.array( [ ] )
	mu1_0 = np.array( [ ] )
	mu2_0 = np.array( [ ] )
	mu1e_0 = np.array( [ ] )
	mu2e_0 = np.array( [ ] )
	a1_0 = np.array( [ ] )
	a2_0 = np.array( [ ] )
	a1e_0 = np.array( [ ] )
	a2e_0 = np.array( [ ] )

	print 'stdev fitting...'

	for index, b in enumerate(b_all):
		print 'bin', index
		N, B = np.histogram(b[ : , 2], bins = 20)
		E = np.sqrt(N)
		E = np.where(E > 0., E, 2.)
		bincenters = 0.5 * (B[1:] + B[:-1])
		opt, cov = curve_fit(func0, bincenters, N, p0 = p0[index], sigma = E, maxfev = 100000)
		
		lowmag = edges[index]
		highmag = edges[index + 1]
		medmag = 0.5 * (lowmag + highmag)

		sigma1, sigma2, m1, m2, a1, a2 = opt[:]
		fiterrors = np.sqrt(np.diagonal(cov))
		sigma1e_0 = np.append(sigma1e_0, fiterrors[0])
		sigma2e_0 = np.append(sigma2e_0, fiterrors[1])
		
		sigma1_0 = np.append(sigma1_0, sigma1)
		sigma2_0 = np.append(sigma2_0, sigma2)
		mag = np.append(mag, medmag)
		print np.column_stack((opt, fiterrors))
		
		x = np.linspace(B[0], B[-1], num = 200)
		test0 = func0(x, opt[0], opt[1], opt[2], opt[3], opt[4], opt[5])
		
		# -----
		# Now plot those histograms, and their constituents
		# -----
		
		plt.plot(x, gauss(x, opt[0], opt[2], opt[4]), linestyle = ':', color = 'k')
		plt.plot(x, gauss(x, opt[1], opt[3], opt[5]), linestyle = '--', color = 'k')
		
		plt.plot(x, test0, color = 'g')
		plt.scatter(bincenters, N, color = 'g', linewidth = 0)
		plt.errorbar(bincenters, N, yerr = E, fmt = None, marker = None)

		plt.xlabel('color (u - r)')
		plt.ylabel('frequency')
		plt.title('iteration 0, %s < r < %s' % (lowmag, highmag))
		plt.grid(True)
		plt.legend(('S1 gaussian', 'S2 gaussian', 'double gaussian'), loc = 'best')
		plt.show()
		

	sigma1_0 = np.absolute(sigma1_0)	
	sigma2_0 = np.absolute(sigma2_0)

	s1start = (0.278, 0.005, -0.075, -19.20, 1.5)
	s1start_b = (0.298, 0.014, -0.067, -19.90, 0.58)

	s2start = (0.175, 0.007, 0.038, -20., 0.94)
	s2start_b = (0.152, 0.008, 0.044, -19.91, 0.94)

	s1fit, s1fit_e = curve_fit(tanhlin, mag, sigma1_0, sigma = sigma1e_0, p0 = s1start, maxfev = 100000)
	print 'blue:'
	print s1fit
	s2fit, s2fit_e = curve_fit(tanhlin, mag, sigma2_0, sigma = sigma2e_0, p0 = s2start, maxfev = 100000)
	print 'red:'
	print s2fit
	
	s1_0 = tanhlin(np.linspace(r_low, r_high, 200), s1fit[0], s1fit[1], s1fit[2], s1fit[3], s1fit[4])
	S1_0 = tanhlin(mag, s1fit[0], s1fit[1], s1fit[2], s1fit[3], s1fit[4])	# specific points

	s2_0 = tanhlin(np.linspace(r_low, r_high, 200), s2fit[0], s2fit[1], s2fit[2], s2fit[3], s2fit[4])
	S2_0 = tanhlin(mag, s2fit[0], s2fit[1], s2fit[2], s2fit[3], s2fit[4])	# specific points
	
	# plot hand-fit s1 vs Baldry and computer-fit
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s1start[0], s1start[1], s1start[2], s1start[3], s1start[4]), color = 'b', linestyle = ':')
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s1start_b[0], s1start_b[1], s1start_b[2], s1start_b[3], s1start_b[4]), color = 'b', linestyle = '--')
	plt.plot(np.linspace(r_low, r_high, 200), s1_0, color = 'b')

	# plot hand-fit s2 vs Baldry and computer-fit
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s2start[0], s2start[1], s2start[2], s2start[3], s2start[4]), color = 'r', linestyle = ':')
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s2start_b[0], s2start_b[1], s2start_b[2], s2start_b[3], s2start_b[4]), color = 'r', linestyle = '--')
	plt.plot(np.linspace(r_low, r_high, 200), s2_0, color = 'r')

	plt.errorbar(mag, sigma1_0, yerr = sigma1e_0, marker = 'd', c = 'b', fmt = 'd', ecolor = 'b')
	plt.errorbar(mag, sigma2_0, yerr = sigma2e_0, marker = 'x', c = 'r', fmt = 'x', ecolor = 'r')
	plt.xlabel('r magnitude')
	plt.ylabel('standard deviation')
	plt.title('color spread of magnitude bins')
	plt.ylim([0.0, 0.9])
	plt.legend(('S1 hand fit', 'S1 Baldry', 'S1 computer fit', 'S2 hand fit', 'S2 Baldry', 'S2 computer fit'), loc = 'best')
	plt.grid(True)
	plt.show()

	# -----
	# Now we start discarding points...
	# the first seven are really quite good, so we take those
	# -----
	
	s1fit1, s1fit1_e = curve_fit(tanhlin, mag[:7], sigma1_0[:7], sigma = sigma1e_0[:7], p0 = s1start, maxfev = 100000)
	s2fit1, s2fit1_e = curve_fit(tanhlin, mag[:7], sigma2_0[:7], sigma = sigma2e_0[:7], p0 = s2start, maxfev = 100000)
	
	s1_1 = tanhlin(np.linspace(r_low, r_high, 200), s1fit1[0], s1fit1[1], s1fit1[2], s1fit1[3], s1fit1[4])
	S1_1 = tanhlin(mag, s1fit1[0], s1fit1[1], s1fit1[2], s1fit1[3], s1fit1[4])	# specific points

	s2_1 = tanhlin(np.linspace(r_low, r_high, 200), s2fit1[0], s2fit1[1], s2fit1[2], s2fit1[3], s2fit1[4])
	S2_1 = tanhlin(mag, s2fit1[0], s2fit1[1], s2fit1[2], s2fit1[3], s2fit1[4])	# specific points

	# plot hand-fit s1 vs Baldry and computer-fit with limited data
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s1start[0], s1start[1], s1start[2], s1start[3], s1start[4]), color = 'b', linestyle = ':')
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s1start_b[0], s1start_b[1], s1start_b[2], s1start_b[3], s1start_b[4]), color = 'b', linestyle = '--')
	plt.plot(np.linspace(r_low, r_high, 200), s1_1, color = 'b')

	# plot hand-fit s2 vs Baldry and computer-fit with limited data
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s2start[0], s2start[1], s2start[2], s2start[3], s2start[4]), color = 'r', linestyle = ':')
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), s2start_b[0], s2start_b[1], s2start_b[2], s2start_b[3], s2start_b[4]), color = 'r', linestyle = '--')
	plt.plot(np.linspace(r_low, r_high, 200), s2_1, color = 'r')

	plt.errorbar(mag[:7], sigma1_0[:7], yerr = sigma1e_0[:7], marker = 'd', fmt = 'x', color = 'b', ecolor = 'b')
	plt.errorbar(mag[:7], sigma2_0[:7], yerr = sigma2e_0[:7], marker = 'x', fmt = 'x', color = 'r', ecolor = 'r')
	plt.xlabel('r magnitude')
	plt.ylabel('standard deviation')
	plt.title('color spread of magnitude bins (first seven points)')
	plt.ylim([0.0, 0.9])
	plt.legend(('blue cloud hand fit', 'blue cloud Baldry', 'blue cloud computer fit', 'red series hand fit', 'red series Baldry', 'red series computer fit'), loc = 'best')
	plt.grid(True)
	plt.show()

	# for the time being, set final sigma equal to the seven-point value

	sigma1_f = np.asfarray(S1_0, dtype=float)
	sigma2_f = np.asfarray(S2_0, dtype=float)
	magnew = mag[:7]
	b_new = [ b0_0, b1_0, b2_0, b3_0, b4_0, b5_0, b6_0 ]

	# define a new function that uses a specified bin
	def func1_0(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[0]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[0]**2.)))
	def func1_1(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[1]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[1]**2.)))
	def func1_2(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[2]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[2]**2.)))
	def func1_3(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[3]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[3]**2.)))
	def func1_4(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[4]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[4]**2.)))
	def func1_5(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[5]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[5]**2.)))
	def func1_6(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[6]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[6]**2.)))
	def func1_7(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[7]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[7]**2.)))
	def func1_8(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[8]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[8]**2.)))
	def func1_9(x, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[9]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[9]**2.)))

	def func1_n(x, n, mu1, mu2, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - mu1)**2. / (2. * sigma1_f[n]**2.))) + a2 * np.exp(-((x - mu2)**2. / (2. * sigma2_f[n]**2.)))

	p1 = np.array( [ [0, 1.8, 2.55, 80., 245.], [1, 2.1, 2.55, 325., 925.], [2, 2.1, 2.55, 325., 925.], [3, 1.8, 2.5, 7500., 20000.], [4, 1.6, 2.55, 15000., 22000.], [5, 1.4, 2.5, 10000., 8000.], [6, 1.5, 2.45, 9000., 10000.] ] )	# new starting values, identical to before, except for elimination of 1st two values
	
	# -----
	# now proceed to second loop, optimizing mean (only for first seven magnitude bins)
	# -----

	print 'mu fitting...'

	for index, b in enumerate(b_new):
		print 'bin', index
		# print b
		N, B = np.histogram(b[ : , 2], bins = 20)
		E = np.sqrt(N)
		E = np.where(E > 0., E, 2.)
		# print 'E:'
		# print E
		bincenters = 0.5 * (B[1:] + B[:-1])
		x = np.linspace(B[0], B[-1], num = 200)
		if index == 0:
			opt, cov = curve_fit(func1_0, bincenters, N, p0 = (p1[0, 1], p1[0, 2], p1[0, 3], p1[0, 4]), sigma = E, maxfev = 100000)
			test1 = func1_0(x, opt[0], opt[1], opt[2], opt[3])
		if index == 1:
			opt, cov = curve_fit(func1_1, bincenters, N, p0 = (p1[1, 1], p1[1, 2], p1[1, 3], p1[1, 4]), sigma = E, maxfev = 100000)
			test1 = func1_1(x, opt[0], opt[1], opt[2], opt[3])
		if index == 2:
			opt, cov = curve_fit(func1_2, bincenters, N, p0 = (p1[2, 1], p1[2, 2], p1[2, 3], p1[2, 4]), sigma = E, maxfev = 100000)
			test1 = func1_2(x, opt[0], opt[1], opt[2], opt[3])
		if index == 3:
			opt, cov = curve_fit(func1_3, bincenters, N, p0 = (p1[3, 1], p1[3, 2], p1[3, 3], p1[3, 4]), sigma = E, maxfev = 100000)
			test1 = func1_3(x, opt[0], opt[1], opt[2], opt[3])
		if index == 4:
			opt, cov = curve_fit(func1_4, bincenters, N, p0 = (p1[4, 1], p1[4, 2], p1[4, 3], p1[4, 4]), sigma = E, maxfev = 100000)
			test1 = func1_4(x, opt[0], opt[1], opt[2], opt[3])
		if index == 5:
			opt, cov = curve_fit(func1_5, bincenters, N, p0 = (p1[5, 1], p1[5, 2], p1[5, 3], p1[5, 4]), sigma = E, maxfev = 100000)
			test1 = func1_5(x, opt[0], opt[1], opt[2], opt[3])
		if index == 6:
			opt, cov = curve_fit(func1_6, bincenters, N, p0 = (p1[6, 1], p1[6, 2], p1[6, 3], p1[6, 4]), sigma = E, maxfev = 100000)
			test1 = func1_6(x, opt[0], opt[1], opt[2], opt[3])

		'''print 'opt: '
		print opt
		print 'cov: '
		print cov'''

		m1, m2, a1, a2 = opt[:]
		# cov is infinite (has to do with passing the index as an argument)
		fiterrors = np.sqrt(np.diagonal(cov))
		mu1e_0 = np.append(mu1e_0, fiterrors[0]) # since mu1 is 2nd diagonal element
		mu2e_0 = np.append(mu2e_0, fiterrors[1]) # since mu2 is 3rd diagonal element
		
		mu1_0 = np.append(mu1_0, m1)
		mu2_0 = np.append(mu2_0, m2)
		

		lowmag = edges[index]
		highmag = edges[index + 1]
		medmag = 0.5 * (lowmag + highmag)
	
		# -----
		# Now plot those histograms, and their constituents
		# -----
		
		'''plt.plot(x, gauss(x, sigma1_f[index], opt[0], opt[2]), linestyle = ':', color = 'k')
		plt.plot(x, gauss(x, sigma2_f[index], opt[1], opt[3]), linestyle = '--', color = 'k')
		
		plt.plot(x, test1, color = 'g')
		plt.scatter(bincenters, N, color = 'g', linewidth = 0)
		plt.errorbar(bincenters, N, yerr = E, fmt = None, marker = None)

		plt.xlabel('color (u - r)')
		plt.ylabel('frequency')
		plt.title('iteration 0, %s < r < %s' % (lowmag, highmag))
		plt.grid(True)
		plt.legend(('blue cloud gaussian', 'red series gaussian', 'double gaussian'), loc = 'best')
		plt.show()'''

	# -----
	# Now fit the means of the gaussians to a tanhlin function
	# -----

	m1start = (1.79, -0.053, -0.363, -20.75, 1.12)
	m1start_b = (1.79, -0.053, -0.363, -20.75, 1.12)

	m2start = (2.279, -0.037, -0.108, -19.81, 0.96)
	m2start_b = (2.279, -0.037, -0.108, -19.81, 0.96)

	mu1fit, mu1fit_e = curve_fit(tanhlin, magnew, mu1_0, p0 = m1start, maxfev = 100000)
	mu2fit, mu2fit_e = curve_fit(tanhlin, magnew, mu2_0, p0 = m2start, maxfev = 100000)

	m1_0 = tanhlin(np.linspace(r_low, r_high, 200), mu1fit[0], mu1fit[1], mu1fit[2], mu1fit[3], mu1fit[4])
	M1_0 = tanhlin(mag, mu1fit[0], mu1fit[1], mu1fit[2], mu1fit[3], mu1fit[4])	# specific points
	m2_0 = tanhlin(np.linspace(r_low, r_high, 200), mu2fit[0], mu2fit[1], mu2fit[2], mu2fit[3], mu2fit[4])
	M2_0 = tanhlin(mag, mu2fit[0], mu2fit[1], mu2fit[2], mu2fit[3], mu2fit[4])	# specific points

	plt.plot(np.linspace(r_low, r_high, 200), m1_0, color = 'b', linestyle = '-')
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), m1start_b[0], m1start_b[1], m1start_b[2], m1start_b[3], m1start_b[4]), color = 'b', linestyle = '--')
	plt.plot(np.linspace(r_low, r_high, 200), m2_0, color = 'r', linestyle = '-')
	plt.plot(np.linspace(r_low, r_high, 200), tanhlin(np.linspace(r_low, r_high, 200), m2start_b[0], m2start_b[1], m2start_b[2], m2start_b[3], m2start_b[4]), color = 'r', linestyle = '--')

	plt.errorbar(magnew, mu1_0, mu1e_0, marker = 'd', fmt = 'd', c = 'b', ecolor = 'b')
	plt.errorbar(magnew, mu2_0, mu2e_0, marker = 'd', fmt = 'd', c = 'r', ecolor = 'r')
	plt.xlabel('r magnitude')
	plt.ylabel('mean color')
	plt.title('color of observed galaxies')
	plt.grid(True)
	plt.ylim([1.0, 2.8])
	plt.legend(('blue cloud computer fit', 'blue cloud Baldry fit', 'red series computer fit', 'red series Baldry fit'), loc = 'best')
	plt.show()

	# define a new function that uses a specified bin
	def func2_0(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[0])**2. / (2. * sigma1_f[0]**2.))) + a2 * np.exp(-((x - M2_0[0])**2. / (2. * sigma2_f[0]**2.)))
	def func2_1(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[1])**2. / (2. * sigma1_f[1]**2.))) + a2 * np.exp(-((x - M2_0[1])**2. / (2. * sigma2_f[1]**2.)))
	def func2_2(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[2])**2. / (2. * sigma1_f[2]**2.))) + a2 * np.exp(-((x - M2_0[2])**2. / (2. * sigma2_f[2]**2.)))
	def func2_3(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[3])**2. / (2. * sigma1_f[3]**2.))) + a2 * np.exp(-((x - M2_0[3])**2. / (2. * sigma2_f[3]**2.)))
	def func2_4(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[4])**2. / (2. * sigma1_f[4]**2.))) + a2 * np.exp(-((x - M2_0[4])**2. / (2. * sigma2_f[4]**2.)))
	def func2_5(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[5])**2. / (2. * sigma1_f[5]**2.))) + a2 * np.exp(-((x - M2_0[5])**2. / (2. * sigma2_f[5]**2.)))
	def func2_6(x, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[6])**2. / (2. * sigma1_f[6]**2.))) + a2 * np.exp(-((x - M2_0[6])**2. / (2. * sigma2_f[6]**2.)))
	
	def func2_n(x, n, a1, a2):
		import numpy as np
		return a1 * np.exp(-((x - M1_0[n])**2. / (2. * sigma1_f[n]**2.))) + a2 * np.exp(-((x - M2_0[n])**2. / (2. * sigma2_f[n]**2.)))

	p2 = np.array( [ [0, 80., 245.], [1, 325., 925.], [2, 325., 925.], [3, 7500., 20000.], [4, 15000., 22000.], [5, 10000., 8000.], [6, 9000., 10000.] ] )	# new starting values, identical to before, except for elimination of 1st two values
	
	# -----
	# now proceed to second loop, optimizing mean (only for first seven magnitude bins)
	# -----

	print 'amplitude fitting...'

	for index, b in enumerate(b_new):
		print 'bin', index
		# print b
		N, B = np.histogram(b[ : , 2], bins = 20)
		E = np.sqrt(N)
		E = np.where(E > 0., E, 2.)
		# print 'E:'
		# print E
		bincenters = 0.5 * (B[1:] + B[:-1])
		x = np.linspace(B[0], B[-1], num = 200)
		if index == 0:
			opt, cov = curve_fit(func2_0, bincenters, N, p0 = (p2[0, 1], p2[0, 2]), sigma = E, maxfev = 100000)
			test1 = func2_0(x, opt[0], opt[1])
		if index == 1:
			opt, cov = curve_fit(func2_1, bincenters, N, p0 = (p2[1, 1], p2[1, 2]), sigma = E, maxfev = 100000)
			test1 = func2_1(x, opt[0], opt[1])
		if index == 2:
			opt, cov = curve_fit(func2_2, bincenters, N, p0 = (p2[2, 1], p2[2, 2]), sigma = E, maxfev = 100000)
			test1 = func2_2(x, opt[0], opt[1])
		if index == 3:
			opt, cov = curve_fit(func2_3, bincenters, N, p0 = (p2[3, 1], p2[3, 2]), sigma = E, maxfev = 100000)
			test1 = func2_3(x, opt[0], opt[1])
		if index == 4:
			opt, cov = curve_fit(func2_4, bincenters, N, p0 = (p2[4, 1], p2[4, 2]), sigma = E, maxfev = 100000)
			test1 = func2_4(x, opt[0], opt[1])
		if index == 5:
			opt, cov = curve_fit(func2_5, bincenters, N, p0 = (p2[5, 1], p2[5, 2]), sigma = E, maxfev = 100000)
			test1 = func2_5(x, opt[0], opt[1])
		if index == 6:
			opt, cov = curve_fit(func2_6, bincenters, N, p0 = (p2[6, 1], p2[6, 2]), sigma = E, maxfev = 100000)
			test1 = func2_6(x, opt[0], opt[1])

		'''print 'opt: '
		print opt
		# print 'cov: '
		# print cov'''

		a1, a2 = opt[:]

		fiterrors = np.sqrt(np.diagonal(cov))
		a1e_0 = np.append(a1e_0, fiterrors[0]) # since mu1 is 2nd diagonal element
		a2e_0 = np.append(a2e_0, fiterrors[1]) # since mu2 is 3rd diagonal element
		
		a1_0 = np.append(a1_0, a1)
		a2_0 = np.append(a2_0, a2)
		
		'''print 'errors:'
		print fiterrors'''

		lowmag = edges[index]
		highmag = edges[index + 1]
		medmag = 0.5 * (lowmag + highmag)

		# -----
		# Now plot those histograms, and their constituents
		# -----
		
		'''plt.plot(x, gauss(x, sigma1_f[index], mu1_0[index], opt[0]), linestyle = ':', color = 'k')
		plt.plot(x, gauss(x, sigma2_f[index], mu2_0[index], opt[1]), linestyle = '--', color = 'k')
		
		plt.plot(x, test1, color = 'g')
		plt.scatter(bincenters, N, color = 'g', linewidth = 0)
		plt.errorbar(bincenters, N, yerr = E, fmt = None, marker = None)

		plt.xlabel('color (u - r)')
		plt.ylabel('frequency')
		plt.title('iteration 0, %s < r < %s' % (lowmag, highmag))
		plt.grid(True)
		plt.legend(('blue cloud gaussian', 'red series gaussian', 'double gaussian'), loc = 'best')
		plt.show()'''

	# -----
	# Now plot the abundance of colors
	# -----

	plt.errorbar(magnew, a1_0, a1e_0, marker = 'd', fmt = 'd', c = 'b', ecolor = 'b')
	plt.errorbar(magnew, a2_0, a2e_0, marker = 'd', fmt = 'd', c = 'r', ecolor = 'r')
	plt.xlabel('r magnitude')
	plt.ylabel('counts')
	plt.title('color abundance of observed galaxies')
	plt.grid(True)
	plt.show()

cmfit('../../../Downloads/Zoo1_CM1_zjpace.fit', 13, 14, 0., 3.5, -24., -15.)
