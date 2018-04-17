from time import time
from pprint import pprint
from dyn_tools.common import *
from dyn_tools.vec import vec
from dyn_tools.gwp import gwp
from dyn_tools.basis import basis
from eig_tools.eig_tools import sinvert
import matplotlib.pyplot as plt


x = gwp( 1, 0, 0, 1)
y = gwp( 1, 10, 0, 1)

b = basis( 0.05, [x] )

start = time()

print( b )
print( )

#b.create_start( )
#print()


N = []
npts = 100
for i in range( npts ):
	b.dynamics( i, 0.001 )
	N.append( b.norm() )

print( b )
print( b.norm() )
print()

plt.plot( range( npts ), [ _.real for _ in N ], 'o' )
plt.show()

print( "CPU time = {0}".format( time() - start ) )
