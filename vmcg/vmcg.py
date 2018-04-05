import numpy as np
from itertools import product
from functools import partial

#________________________________________________________________
# Potential energy
#________________________________________________________________
def potential( r ):
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
		if ( type(other) == type(self) ):
			return vec( self.x + other.x, self.y + other.y, self.z + other.z )
		else: return vec( self.x + other, self.y + other, self.z + other )

	def __radd__( self, other ):
		return vec( other + self.x, other + self.y, other + self.z)

	def __sub__ ( self, other ):
		if ( type(other) == type(self) ):
			return vec( self.x - other.x, self.y - other.y, self.z - other.z )
		else: return vec( self.x - other, self.y - other, self.z - other )

	def __rsub__( self, other ):
		return vec( other - self.x, other - self.y, other - self.z)

	def scalar( self, other ):
		return self.x * other.x + self.y * other.y + self.z * other.z 

	def __mul__ ( self, other ):
		if ( type(other) == type(self) ):
			return vec( self.x * other.x, self.y * other.y, self.z * other.z )
		else: return vec( other * self.x, other * self.y, other * self.z )

	def __rmul__( self, other ):
		return vec( other * self.x, other * self.y, other * self.z )

	def __pow__ ( self, other ):
		return vec( self.x**other, self.y**other, self.z**other )

	def __truediv__( self, other ):
		return vec( self.x / other, self.y / other, self.z / other )

	def call( self, i ):
		if i == 0: return self.x
		if i == 1: return self.y
		if i == 2: return self.z

	def sum( self ):
		return self.x + self.y + self.z	

	def conj( self ):
		return vec( np.conj( self.x ), np.conj( self.y ), np.conj( self.z ) )

	def exp( self ):
		return np.exp( self.x + self.y + self.z )  

	def min( self, other ):
		return vec( min( self.x, other.x ), min( self.y, other.y ), min( self.z, other.z ) )

	def max( self, other ):
		return vec( max( self.x, other.x ), max( self.y, other.y ), max( self.z, other.z ) )

	def __str__( self ):
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
		return np.conj( self.D ) * other.D * ( np.pi / self.omega )**( 1.5 ) * \
		       ( 0.5 * ( self.xi.conj() + other.xi )**2 / \
			       ( self.omega + other.omega ) + \
			         self.eta.conj() + other.eta ).exp()

	def __rmul__ ( self, other ):
		return GWP( other * self.D, self.q, self.p, self.omega )

	def der1( self, other ):
		return   self * other * ( self.xi.conj() + other.xi ) / ( self.omega + other.omega )

	def der11( self, other ):
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
		return  self * other * ( 1.0 / ( self.omega + other.omega ) + \
					 ( self.xi.conj() + other.xi )**2 / \
					 ( self.omega + other.omega )**2 )

	def psi( self, r ):
		return np.exp( -self.omega * r.scalar( r ) + self.xi.scalar( r ) + self.eta.sum() )

	def T( self, other ):
		return  0.5 * self * other * ( 3.0 * other.omega - other.xi.scalar( other.xi ) ) + \
			other.omega * other.xi.scalar( self.der1( other ) ) - \
			0.5 * self.omega**2 * self.der2( other ).sum()

	def V( self, other, NPTS ):
		x, w = 	np.polynomial.legendre.leggauss( NPTS )
		v = 0.0
		lb = (self.q - 5.0 / self.omega).min( other.q - 5.0 / self.omega )
		rb = (self.q + 5.0 / self.omega).max( other.q + 5.0 / self.omega )
		shift_x = partial( lambda a, b, x: 0.5 * ( b - a ) * x + 0.5 * ( a + b ), \
				   lb, rb )

		x_corr = list( map( shift_x, x ) )
		
		for i, j, k in product( range( NPTS ), range( NPTS ), range( NPTS ) ):
			r = vec( x_corr[i], x_corr[j], x_corr[k] ) 
			v += w[i] * w[j] * w[k] * np.conj( self.psi( r ) ) * potential( r ) * other.psi( r )
		return v

	def __str__( self ):
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

print( x.V( x, 48 ) )


Nfunc = 1
Ndim = 1
omega = 1
D = [ 1.0 ]
q = [ 1.0 ]
p = [ 0.0 ]

S = np.zeros( ( Nfunc, Nfunc ) )
S_alpha0 = np.zeros( ( Ndim * Nfunc, Nfunc ) )
S_0alpha = np.zeros( ( Nfunc, Ndim * Nfunc ) )
S_alphabeta = np.zeros( ( Ndim * Nfunc, Ndim * Nfunc ) )
H = np.zeros( ( Nfunc, Nfunc ) )
H_alpha0 = np.zeros( ( Ndim * Nfunc, Nfunc ) )


