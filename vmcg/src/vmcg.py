from time import time
from pprint import pprint
from dyn_tools.common import *
from dyn_tools.vec import vec
from dyn_tools.gwp import gwp
from dyn_tools.basis import basis
from eig_tools.eig_tools import sinvert
import matplotlib.pyplot as plt


x = gwp( 1, 1, 0, 1)
y = gwp( 1, 10, 0, 1)

b = basis( 0.05, [x] )

start = time()

print( b )
print( b.norm() )
print( b.energy() )
print( )

N = []
E = []
q = []
p = []
npts = 100
for i in range( npts ):
	b.dynamics( 0.001 )
	N.append( b.norm() )
	E.append( b.energy() )
	q.append( b.basis[0].q.call(0) )
	p.append( b.basis[0].p.call(0) )

print( b )
print( b.norm() )
print( b.energy() )
print()

#plt.plot( range( npts ), [ _.real for _ in E ], 'o' )
#plt.plot( q, p )
#plt.show()

print( "CPU time = {0}".format( time() - start ) )
