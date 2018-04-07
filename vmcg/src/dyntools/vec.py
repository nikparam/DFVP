from .common import * 
#________________________________________________________________
# 3d vector
#________________________________________________________________

class vec():
	def __init__ ( self, coords = None ):
		if coords is None:
			self.coords = [ 0.0 ]
		elif isinstance( coords, vec ):
			self.coords = coords.coords[:]
		elif type( coords ) != type( [] ):
			self.coords = [ coords ]
		elif not coords:
			self.coords = [ 0.0 ]
		else:
			self.coords = coords[:]
			
	def __add__ ( self, other ):
		'''
		  vec(a) + vec(b) = vec( a.x + b.x, a.y + b.y, a.z + b.z )
		'''

		if ( type(other) == type(self) ):
			return vec( lmap( lambda x, y: x + y, self.coords, other.coords )  )
		else: return vec( lmap( lambda x: x + other, self.coords ) )

	def __iadd__ ( self, other ):

		return self + other

	def __radd__( self, other ):
		'''
		  num + vec(a) = vec( num + a.x, num + a.y, num + a.z )
		'''

		return vec( lmap( lambda x: x + other, self.coords ) )

	def __sub__ ( self, other ):
		'''
		if we subtract one vec from another:
		  vec(a) - vec(b) = vec( a.x - b.x, a.y - b.y, a.z - b.z )
		if we subtract number from vector:
		  vec(a) - num = vec( a.x - num, a.y - num, a.z - num )
		'''

		if ( type(other) == type(self) ):
			return vec( lmap( lambda x, y: x - y, self.coords, other.coords ) )
		else: return vec( lmap( lambda x: x - other, self.coords ) )

	def __rsub__( self, other ):
		'''
		  num - vec(a) = vec( num - a.x, num - a.y, num - a.z )
		'''

		return vec( lmap( lambda x: other - x, self.coords )  )

	def __neg__( self ):
		return vec( lmap( lambda x: -x, self.coords ) )

	def sum( self ):
		'''
		sum of elements
		'''

		return reduce( lambda a, x: a + x, self.coords )	

	def scalar( self, other ):
		'''
		dot product of two vectors
		'''

		return ( self * other ).sum() 

	def __mul__ ( self, other ):
		'''
		if we multiply two vectors:
		 vec(a) * vec(b) = vec( a.x * b.x, a.y * b.y, a.z * b.z ) Hadamard multiplication
		if we multiply vector by number:
		 vec(a) * num = vec( a.x * num, a.y * num, a.z * num )
		'''

		if ( type(other) == type(self) ):
			return vec( lmap( lambda x, y: x * y, self.coords, other.coords ) )
		else: return vec( lmap( lambda x: x * other, self.coords ) )

	def __rmul__( self, other ):
		'''
		 vec(a) * num = vec( a.x * num, a.y * num, a.z * num )
		'''

		return vec( lmap( lambda x: x * other, self.coords ) )

	def __pow__ ( self, other ):
		'''
		 vec(a) ** num = vec( a.x ** num, a.y ** num, a.z ** num )
		'''
		
		return vec( lmap( lambda x: x ** other, self.coords ) )

	def __rtruediv__( self, other ):

		return vec( lmap( lambda x: other / x, self.coords ) )

	def __truediv__( self, other ):
		'''
		 vec(a) / num = vec( a.x / num, a.y / num, a.z / num )
		'''
		if ( type( self ) != type( other ) ): 
			return vec( lmap( lambda x: x / other, self.coords ) )
		else:
			return vec( lmap( lambda x, y: x / y, self.coords, other.coords ) )

	def __lt__( self, other ):

		if True in lmap( lambda x, y: x < y, self.coords, other.coords ):
			return True
		else:
			return False

	def __le__( self, other ):

		if False not in lmap( lambda x, y: x <= y, self.coords, other.coords ):
			return True
		else:
			return False

	def __ge__( self, other ):

		if ( self <= other ):
			return False
		else:
			return True

	def __gt__( self, other ):

		if ( self < other ):
			return False
		else:
			return True

	def __eq__( self, other ):

		if False not in lmap( lambda x, y: x == y, self.coords, other.coords ):
			return True
		else: 
			return False

	def __ne__( self, other ):

		if ( self == other ):
			return False
		else: 
			return True

	def mult( self ):
		return reduce( lambda a, x: a * x, self.coords, 1 )

	def conj( self ):
		'''
		complex conjugate of vector
		'''

		return vec( lmap( lambda x: np.conj( x ), self.coords ) )

	def exp( self ):
		'''
		exp( vec(a) ) = vec( exp( vec.coords ) )
		'''

		return vec( lmap( lambda x: np.exp( x ), self.coords ) )

	def log( self ):
		'''
		exp( vec(a) ) = vec( exp( vec.coords ) )
		'''

		return vec( lmap( lambda x: np.log( x ), self.coords ) )

	def min( self, other ):
		'''
		compare elements of two vecs, build vec, that consists of min coords
		'''

		return vec( lmap( lambda x, y: min(x,y), self.coords, other.coords ) )

	def max( self, other ):
		'''
		compare elements of two vecs, build vec, that consists of max coords
		'''

		return vec( lmap( lambda x, y: max(x,y), self.coords, other.coords ) )

	def size( self ):

		return len( self.coords )

	def call( self, i ):

		return self.coords[i]

	def __str__( self ):
		'''
		print vec
		'''

		return reduce( lambda a, x: a + str( x ) + ' ', self.coords, '' )

