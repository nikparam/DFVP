import numpy as np
from itertools import product
from functools import partial, reduce

def compose( *func ):
	init, *rest = reversed(func)
	return lambda *args, **kws: reduce( lambda a, x: x(a), rest, init(*args, **kws) )

lmap = compose( list, map )
mtrx = partial( np.zeros, dtype = complex ) 
hermite = compose( np.conj, np.transpose )

def matrix( func ):
	return lambda right, left: \
				   lmap(\
					 lambda f: lmap( f, right ), \
					 lmap(\
					       lambda x: partial( func, x ), \
					       left \
					     )\
				        )

def print_matrix( matrix ):
	string = ''
	for row in matrix:
		if type( row ) == type( matrix ):
			for _ in row:
				string += '({:.4f})\t'.format( _ )
			string += '\n'
		else:
			string += '({:.4f})\n'.format( row )
	return string

