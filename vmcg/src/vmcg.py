from time import time
from pprint import pprint
from dyntools.common import *
from dyntools.vec import vec
from dyntools.gwp import gwp
from dyntools.basis import basis
from eig_tools.eig_tools import sinvert

x = gwp( 1, 0, 0, 1)
y = gwp( 1, 1, 1, 1)

b = basis( 0.05 )
print( b )
print( )

b.create_start()
print( print_matrix( b.xmatrix() ) ) 
