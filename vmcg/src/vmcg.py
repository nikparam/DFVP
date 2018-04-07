from time import time
from pprint import pprint
from dyn_tools.common import *
from dyn_tools.vec import vec
from dyn_tools.gwp import gwp
from dyn_tools.basis import basis
from eig_tools.eig_tools import sinvert

x = gwp( 1, 0, 0, 1)
y = gwp( 1, 1, 0, 1)

b = basis( 0.05, [x] )
print( b )
print( )

b.create_start()
b.dynamics( 0.01 )
print( b )

