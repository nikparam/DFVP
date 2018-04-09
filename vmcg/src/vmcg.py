from time import time
from pprint import pprint
from dyn_tools.common import *
from dyn_tools.vec import vec
from dyn_tools.gwp import gwp
from dyn_tools.basis import basis
from eig_tools.eig_tools import sinvert

x = gwp( 1, 0, 0, 1)
y = gwp( 1, 10, 0, 1)

b = basis( 0.05, [x] )

start = time()

print( b )
print( )

b.create_start( )
print()
for _ in range(1000):
	b.dynamics( 1e-10 )

print( b )
print( b.norm() )
print( )

print( "CPU time = {0}".format( time() - start ) )
