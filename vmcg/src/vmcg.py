import numpy as np
from itertools import product
from functools import partial
from time import time
from mean_v import mean_v 



#________________________________________________________________
# Potential energy
#________________________________________________________________
def potential( r ):
	'''
	potential energy surface 
	for <g_i|V|g_j> calculation
	'''

	return 0.5 * ( r.x**2 + r.y**2 + r.z**2 )

#________________________________________________________________
# 3d vector
#________________________________________________________________

class vec():
	def __init__ ( self, x, y, z ):
		self.x = x
		self.y = y
		self.z = z

	def __add__ ( self, other ):
		'''
		  vec(a) + vec(b) = vec( a.x + b.x, a.y + b.y, a.z + b.z )
		'''

		if ( type(other) == type(self) ):
			return vec( self.x + other.x, self.y + other.y, self.z + other.z )
		else: return vec( self.x + other, self.y + other, self.z + other )

	def __radd__( self, other ):
		'''
		  num + vec(a) = vec( num + a.x, num + a.y, num + a.z )
		'''

		return vec( other + self.x, other + self.y, other + self.z)

	def __sub__ ( self, other ):
		'''
		if we subtract one vec from another:
		  vec(a) - vec(b) = vec( a.x - b.x, a.y - b.y, a.z - b.z )
		if we subtract number from vector:
		  vec(a) - num = vec( a.x - num, a.y - num, a.z - num )
		'''

		if ( type(other) == type(self) ):
			return vec( self.x - other.x, self.y - other.y, self.z - other.z )
		else: return vec( self.x - other, self.y - other, self.z - other )

	def __rsub__( self, other ):
		'''
		  num - vec(a) = vec( num - a.x, num - a.y, num - a.z )
		'''

		return vec( other - self.x, other - self.y, other - self.z)

	def scalar( self, other ):
		'''
		dot product of two vectors
		'''

		return self.x * other.x + self.y * other.y + self.z * other.z 

	def __mul__ ( self, other ):
		'''
		if we multiply two vectors:
		 vec(a) * vec(b) = vec( a.x * b.x, a.y * b.y, a.z * b.z ) Hadamard multiplication
		if we multiply vector by number:
		 vec(a) * num = vec( a.x * num, a.y * num, a.z * num )
		'''

		if ( type(other) == type(self) ):
			return vec( self.x * other.x, self.y * other.y, self.z * other.z )
		else: return vec( other * self.x, other * self.y, other * self.z )

	def __rmul__( self, other ):
		'''
		 vec(a) * num = vec( a.x * num, a.y * num, a.z * num )
		'''

		return vec( other * self.x, other * self.y, other * self.z )

	def __pow__ ( self, other ):
		'''
		 vec(a) ** num = vec( a.x ** num, a.y ** num, a.z ** num )
		'''
		
		return vec( self.x**other, self.y**other, self.z**other )

	def __truediv__( self, other ):
		'''
		 vec(a) / num = vec( a.x / num, a.y / num, a.z / num )
		'''
		
		return vec( self.x / other, self.y / other, self.z / other )

	def call( self, i ):
		'''
		cheap way to call element of vec
		'''

		if i == 0: return self.x
		if i == 1: return self.y
		if i == 2: return self.z

	def sum( self ):
		'''
		sum of elements
		'''

		return self.x + self.y + self.z	

	def conj( self ):
		'''
		complex conjugate of vector
		'''

		return vec( np.conj( self.x ), np.conj( self.y ), np.conj( self.z ) )

	def exp( self ):
		'''
		exp( vec(a) ) = exp( sum( vec ) )
		'''

		return np.exp( self.sum() )  

	def min( self, other ):
		'''
		compare elements of two vecs, build vec, that consists of min coords
		'''

		return vec( min( self.x, other.x ), min( self.y, other.y ), min( self.z, other.z ) )

	def max( self, other ):
		'''
		compare elements of two vecs, build vec, that consists of max coords
		'''

		return vec( max( self.x, other.x ), max( self.y, other.y ), max( self.z, other.z ) )

	def list( self ):
		'''
		turn vec to list
		'''
		return [ self.x, self.y, self.z ]

	def __str__( self ):
		'''
		print vec
		'''

		return "[ {0}, {1}, {2} ]".format( self.x, self.y, self.z )

#________________________________________________________________
# Gaussian wave packet
#________________________________________________________________

class GWP( ):
	def __init__( self, coeff = 1.0, q = None, p = None, omega = 1.0 ):
		self.D = coeff

		if q is None or q == []:
			self.q = vec( 0.0, 0.0, 0.0 )
		elif isinstance( q, vec ) == True:
			self.q = vec( q.x, q.y, q.z )
		elif type( q ) != type( [] ):
			self.q = vec( q, 0.0, 0.0 )
		else:
			self.q = vec( q[0], q[1], q[2] )

		if p is None or p == []:
			self.p = vec( 0.0, 0.0, 0.0 )
		elif isinstance( p, vec ) == True:
			self.p = vec( p.x, p.y, p.z )
		elif type( p ) != type( [] ):
			self.p = vec( p, 0.0, 0.0 )
		else:
			self.p = vec( p[0], p[1], p[2] )
	
		self.omega = omega
		self.xi = self.omega * self.q + 1j * self.p
		self.eta = 0.25 * ( np.log( self.omega / np.pi ) - 2 * self.omega * self.q**2 ) - 1j * self.q * self.p

	def __mul__( self, other ):
		'''
		overlap <g_i|g_j>
		'''

		return np.conj( self.D ) * other.D * ( np.pi / self.omega )**( 1.5 ) * \
		       ( 0.5 * ( self.xi.conj() + other.xi )**2 / \
			       ( self.omega + other.omega ) + \
			         self.eta.conj() + other.eta ).exp()

	def __rmul__ ( self, other ):
		'''
		create GWP with coeff D_1 = num * D_0
		'''

		return GWP( other * self.D, self.q, self.p, self.omega )

	def der1( self, other ):
		'''
		      d  g_j 
		<g_i| ------ >
		      d xi_j
		'''

		return   self * other * ( self.xi.conj() + other.xi ) / ( self.omega + other.omega )

	def der11( self, other ):
		'''
		  d  g_i | d  g_j 
		< ------ | ------ >
		  d xi_i | d xi_j
		'''

		s = np.zeros( (3,3), dtype = complex )
		N = self * other
		for i in range(3):
			s[i][i] = N * ( 1.0 / ( self.omega + other.omega ) + \
					( self.xi.conj().call(i) + other.xi.call(i) )**2  / ( self.omega + other.omega )**2 )
			for j in range( i+1, 3 ):
				s[i][j] = N * ( self.xi.conj().call(i) + other.xi.call(i) ) * \
					      ( self.xi.conj().call(j) + other.xi.call(j) ) / ( self.omega + other.omega )**2
		return s		
		
		

	def der2( self, other ):
		'''
		      d^2  g_j 
		<g_i| -------- >
		      d xi_j^2
		'''

		return  self * other * ( 1.0 / ( self.omega + other.omega ) + \
					 ( self.xi.conj() + other.xi )**2 / \
					 ( self.omega + other.omega )**2 )

	def psi( self, r ):
		'''
		psi = exp( -0.5 * omega * r^2 + xi * r + eta )
		'''

		return ( -0.5 * self.omega * r * r + self.xi * r + self.eta ).exp()

	def T( self, other ):
		'''
			                     d^2    d^2    d^2
		<g_i|T|g_j> = -0.5 * < g_i | ---- + ---- + ---- | g_j >
			                     dx^2   dy^2   dz^2
		'''

		return  0.5 * self * other * ( 3.0 * other.omega - other.xi.scalar( other.xi ) ) + \
			other.omega * other.xi.scalar( self.der1( other ) ) - \
			0.5 * self.omega**2 * self.der2( other ).sum()


	def V( self, other, NPTS ):
		'''
		<g_i|V|g_j> calculated with Gauss Legendre quadratures
		'''

		x, w = 	np.polynomial.legendre.leggauss( NPTS )
		v = 0.0
		lb = (self.q - 6.0 / self.omega).min( other.q - 6.0 / self.omega )
		rb = (self.q + 6.0 / self.omega).max( other.q + 6.0 / self.omega )

		shift = lambda x, a, b: 0.5 * ( b - a ) * x + 0.5 * ( a + b )

		for i, j, k in product( range( NPTS ), range( NPTS ), range( NPTS ) ):
			r = vec( shift( x[i], lb.x, rb.x ), \
				 shift( x[j], lb.y, rb.y ), \
				 shift( x[k], lb.z, rb.z ) ) 
			v += 0.125 * w[i] * w[j] * w[k] * \
			     ( rb.x - lb.x ) * ( rb.y - lb.y ) * ( rb.z - lb.z ) * \
			     np.conj( self.psi( r ) ) * potential( r ) * other.psi( r )
		return v

	def H( self, other, NPTS ):
		'''
		<g_i|H|g_j> = <g_i|T|g_j> + <g_i|V|g_j>
		<g_i|V|g_j> is calculated with Fortran90 subroutine
		with Gauss Legendre quadratures
		see file mean_v.f90
		'''

		t = self.T( other )
		v = mean_v( NPTS, self.xi.list(), self.eta.list(), self.omega, \
				     other.xi.list(), other.eta.list(), other.omega )
		return t + v

	def __str__( self ):
		'''
		print GWP
		'''

		return  "D = {0:g}".format( self.D ) + ", omega = {0:g}".format( self.omega ) + '\n' + \
			"q = " + str( self.q ) + '\n' + \
			"p = " +  str( self.p )
		
		

class basis():
	def __init__( self, basis = None ):
		if basis is None: self.basis = []
		else: self.basis = basis[:]

	def read_basis( self ):
		try:
			with open('basis.inp','r') as fin:
				N = fin.readline()
				for _ in range( N ):
					string = fin.readline()
					omega = string.split[0]
					coeff = string.split[0]
					q = string.split[1]
					p = string.split[2]
					g = GWP( coeff, q, p, omega )
					self.basis.append( g )

		except:
			print( "Something went wrong" )


x = GWP( 1, 0, 0, 1 )
y = GWP( 1, 1, 1, 1 )

start = time()
print( x.H( x, 50 ) )
print( "fortran CPU time = {0:g}".format( time() - start ) )

start = time()
print( x.T( x ) + x.V( x, 50 ) )
print( "python CPU time = {0:g}".format( time() - start ) )

