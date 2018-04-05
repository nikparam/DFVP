class vec():
	def __init__( self, x, y, z ):
		self.x = x
		self.y = y
		self.z = z

	def __add__( self, other ):
		return vec( self.x + other.x, self.y + other.y, self.z + other.z)

	def __mul__ ( self, other ):
		return self.x * other.x + self.y * other.y + self.z * other.z 

	def __rmul__ ( self, other ):
		return vec( other * self.x, other * self.y, other * self.z )

	def __str__ ( self ):
		return "{0:g}i + {1:g}j + {2:g}k".format( self.x, self.y, self.z )

	def div( self, other ):
		return vec( self.x / other.x, self.y / other.y, self.z / other.z, )

class test( ):
	def __init__( self, m = None ):
		self.m = []
		if m is None or len(m) == 0 : 
			self.m.append( vec( 0.0, 0.0, 0.0 ) )
		else: 
			for _ in m:
				self.m.append( vec( _[0], _[1], _[2] ) )

	def func( self, x, y, z ):
		self.m.append( super().__init__( x, y, z ) )



v1 = vec( 1, 2, 3 )
v2 = vec( 1, 2, 3 )
print( type(v1) == type(v2) )

x = [ [ 1, 2, 3 ], [4, 5, 6] ]
b = test( x )
for _ in b.m:
	print( _ )
