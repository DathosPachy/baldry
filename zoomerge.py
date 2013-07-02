# merge tables for Zoo1 and Zoo2

def loadfiles(z1file, z2file, morph_col):
	import numpy as np
	from pyfits import open
	zoo1_raw = open(z1file)
	zoo1 = zoo1_raw[1].data
	zoo2_raw = open(z2file)
	zoo2 = zoo2_raw[1].data
	zoo1_objID = zoo1.field(0)
	zoo2_objID = zoo2.field(232)
	
	cols1 = zoo1_raw[1].columns
	cols2 = zoo2_raw[1].columns
	cols1names = cols1.names
	cols2names = cols2.names
	
	print 'Building Zoo2 CMD based on column', morph_col, ':', cols2names[morph_col]

	u, r = np.asfarray(zoo1.field(7), dtype = float), np.asfarray(zoo1.field(8), dtype = float)
	Z1 = np.column_stack( ( zoo1_objID, u, r, u - r ) )
	print 'Z1:'
	print Z1

	morph = np.asfarray(zoo2.field(morph_col), dtype = float)
	Z2 = np.column_stack( ( zoo2_objID, morph ) )
	print 'Z2:'
	print Z2
	
	sort_idx = np.argsort(Z1[:, 0])
	print 'sort_idx:'
	print sort_idx
	sorted_insert = np.searchsorted(Z1[sort_idx, 0], Z2[:, 0], side='left')
	print 'sorted_insert:'
	print sorted_insert
	unsorted_insert = np.take(sort_idx, sorted_insert)
	print 'unsorted_insert:'
	print unsorted_insert
	
	zoo2_morph = np.hstack( ( Z2, Z1[unsorted_insert, 1:] ) )
	print zoo2_morph
	return zoo2_morph

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
