# merge tables for Zoo1 and Zoo2

def loadfiles(z1file, z2file, morph_col):
	import numpy as np
	from pyfits import open
	zoo1 = open(z1file)[1].data
	zoo2 = open(z2file)[1].data
	zoo1_objID = zoo1.field(0)
	zoo2_objID = zoo2.field(232)
	
	# zoo1 = zoo1[(zoo1_objID != -9999)]
	# zoo2 = zoo2[(zoo2_objID != -9999)]

	u, r = np.asfarray(zoo1.field(7), dtype = float), np.asfarray(zoo1.field(8), dtype = float)
	print u, r, zoo1_objID
	Z1 = np.hstack( ( zoo1_objID, u, r ) )
	print type(Z1)

	morph = np.asfarray(zoo2.field(morph_col), dtype = float)
	Z2 = np.hstack( ( zoo2_objID, morph ) )
	print Z2
	
	'''sort_idx = np.argsort(Z1[:, 0])
	sorted_insert = np.searchsorted(Z1[:, 0], sorter = sort_idx)
	unsorted_insert_np.take(sort_idx, sorted_insert)
	
	zoo2_morph = np.hstack((Z2, Z1[unsorted_insert, 1:]))'''

'''
def usefits(filename, ucol, rcol, votecol, (cl, ch), (rl, rh)):
	import numpy as np
	from pyfits import open
	hdulist = open(filename)
	raw = hdulist[1].data
	u, r = raw.field(ucol), raw.field(rcol)
	
	colcond = ((u - r) < ch) & ((u - r) > cl) & (r > rl) & (r < rh)
	raw = raw[colcond]
	print raw
	u, r, vote = np.asfarray(raw.field(ucol), dtype = float), np.asfarray(raw.field(rcol), dtype = float), np.asfarray(raw.field(votecol), dtype = float)
	mixed_data = np.column_stack((u, r, u - r, vote)) # use for vote morphology
	mixed_data = np.column_stack((u, r, u - r)) # use for base case
	print 'mixed data:', '\n', mixed_data
	return mixed_data
'''

loadfiles('../../Downloads/Zoo1_CM3_zjpace.fit', '/data/lucifer1.1/users/zpace/gz2_debiased_main.fits', 56)
