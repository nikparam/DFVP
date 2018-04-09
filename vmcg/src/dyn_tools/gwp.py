from .common import *
from mean_v.mean_v import mean_v 
from .vec import vec

#________________________________________________________________
# Potential energy
#________________________________________________________________
def potential( r ):
	'''
	potential energy surface 
	for <g_i|V|g_j> calculation
	'''

#	return 0.0
	return 0.5 * r **2

#________________________________________________________________
# Gaussian wave packet
#________________________________________________________________

class gwp( ):
	def __init__( self, coeff = 1.0, q = None, p = None, omega = 1.0 ):
		self.D = coeff
		self.q = vec( q )	
		self.p = vec( p )		
		self.omega = vec( omega )
		self.xi = self.omega * self.q + 1j * self.p
		self.eta = 0.25 * ( ( self.omega / np.pi ).log() - 2 * self.omega * self.q**2 ) - 1j * self.q * self.p

	def __mul__( self, other ):
		'''
		overlap <g_i|g_j>
		'''

		return ( np.pi / self.omega )**( 0.5 ) * \
		       ( 0.5 * ( self.xi.conj() + other.xi )**2 / \
			       ( self.omega + other.omega ) + \
			         self.eta.conj() + other.eta ).exp()

	def __rmul__ ( self, other ):
		'''
		create GWP with coeff D_1 = num * D_0
		'''

		return gwp( other * self.D, self.q, self.p, self.omega )

	def der1( self, other ):
		'''
		      d  g_j 
		<g_i| ------ >
		      d xi_j
		'''

		return   (  self * other ).mult()* \
			 ( self.xi.conj() + other.xi ) / \
			 ( self.omega + other.omega )

	def der11( self, other ):
		'''
		  d  g_i | d  g_j 
		< ------ | ------ >
		  d xi_i | d xi_j
		'''

		s = mtrx( ( self.xi.size(), other.xi.size() ) )
		N = ( self * other ).mult()
		for i in range( self.xi.size() ):
			s[i][i] = N * ( 1.0 / ( self.omega.call(i) + other.omega.call(i) ) + \
				        ( self.xi.conj().call(i) + other.xi.call(i) )**2  / \
				        ( self.omega.call(i) + other.omega.call(i) )**2 ) 
			for j in range( i+1, other.xi.size() ):
				s[i][j] = N * ( self.xi.conj().call(i) + other.xi.call(i) ) / \
					      ( self.omega.call(i) + other.omega.call(i) ) * \
					      ( self.xi.conj().call(j) + other.xi.call(j) ) / \
					      ( self.omega.call(j) + other.omega.call(j) )
		return s		
		
		

	def der2( self, other ):
		'''
		      d^2  g_j 
		<g_i| -------- >
		      d xi_j^2
		'''

		return  ( self * other ).mult() * \
			( 1.0 / ( self.omega + other.omega ) + \
			( self.xi.conj() + other.xi )**2 / \
			( self.omega + other.omega )**2 )

	def psi( self, r, i ):
		'''
		psi = exp( -0.5 * omega * r^2 + xi * r + eta )
		'''

		return np.exp( -0.5 * self.omega.call(i) * r * r + self.xi.call(i) * r + self.eta.call(i) )

	def T( self, other ):
		'''
			                     d^2    d^2    d^2
		<g_i|T|g_j> = -0.5 * < g_i | ---- + ---- + ---- | g_j >
			                     dx^2   dy^2   dz^2
		'''

		return  0.5 * ( self * other ).mult() * ( other.omega - other.xi * other.xi ) + \
			other.omega * other.xi * self.der1( other ) - \
			0.5 * self.omega**2 * self.der2( other )


	def Vp( self, other, NPTS ):
		'''
		<g_i|V|g_j> calculated with Gauss Legendre quadratures
		'''

		x, w = 	np.polynomial.legendre.leggauss( NPTS )
		lb = ( self.q - 5.0 / self.omega ).min( other.q - 5.0 / self.omega )
		rb = ( self.q + 5.0 / self.omega ).max( other.q + 5.0 / self.omega )

		shift = lambda x, a, b: 0.5 * ( b - a ) * x + 0.5 * ( a + b )

		v = np.zeros( self.q.size(), dtype = complex )

		for i in range( self.q.size() ):
			for j in range( NPTS ):
				r = shift( x[j], lb.call(i), rb.call(i) )
				v[i] += w[j] * np.conj( self.psi( r, i ) ) * potential( r ) * other.psi( r, i )
			v[i] *= 0.5 * ( rb.call(i) - lb.call(i) ) 
		return vec( v )

	def Vf( self, other, NPTS ):
		'''
		<g_i|V|g_j> is calculated with Fortran90 subroutine
		with Gauss Legendre quadratures
		see file mean_v.f90
		'''

		v = mean_v( NPTS, self.xi.coords, self.eta.coords, self.omega.coords, \
				  other.xi.coords, other.eta.coords, other.omega.coords )
		return vec( v )

	def H1( self, other ):

		return ( self.xi.conj() + other.xi ) / \
		       ( self.omega + other.omega ) * self.T( other ) + self.Vp( other, 50 )

	def __str__( self ):
		'''
		print GWP
		'''

		return  "D = {0:g}".format( self.D ) + \
			", omega = " + str( self.omega ) + '\n' + \
			"q = " + str( self.q ) + '\n' + \
			"p = " +  str( self.p )


