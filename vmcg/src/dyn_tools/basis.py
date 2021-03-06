from .common import *
from .vec import vec
from .gwp import gwp
from eig_tools.eig_tools import cholesky, svd, lowdin, sinvert

#________________________________________________________________
# Orthonormal basis
#________________________________________________________________

class basis():
	def __init__( self, E, basis = None, omega = 1.0, Ndim = 1 ):

		if basis is None:
			self.basis = []
			self.basis_tmp = []
			self.Ndim = Ndim
	
			a = 0.2
			b = 0.2
	
			q_max = vec( np.sqrt( 2 * E / omega ) )
			q = vec( [ 0.0 ] * Ndim )
	
			p_max = vec( np.sqrt( 2 * E ) )
			p = vec( [ 0.0 ] * Ndim )
	
			while ( p <= p_max ):
				if  p == vec( [ 0.0 ] * Ndim ) : pass
				else:
					self.basis.append( gwp( 1., vec( [ 0.0 ] * Ndim ), -p, vec( [ omega ] * Ndim ) ) )
					self.basis.append( gwp( 1., vec( [ 0.0 ] * Ndim ), p, vec( [ omega ] * Ndim ) ) )
				p += b
	
			while ( q <= q_max ):
				if q == vec( [ 0.0 ] * Ndim ) : pass
				else:
					self.basis.append( gwp( 1., -q, vec( [ 0.0 ] * Ndim ), vec( [ omega ] * Ndim ) ) )
					self.basis.append( gwp( 1., q, vec( [ 0.0 ] * Ndim ), vec( [ omega ] * Ndim ) ) )
				q += a
		else:
			self.basis = basis
			self.basis_tmp = []
			self.Ndim = Ndim

	def size( self ):
		return len( self.basis )

	def call( self, i ):
		return self.basis[i]

	def smatrix( self ):
		S = mtrx( ( self.size(), self.size() ) )
		for i, j in product( range( self.size() ), range( self.size() ) ):
			S[i][j] = ( self.basis[i] * self.basis[j] ).mult()
		return S

	def tmatrix( self ):
		T = mtrx( ( self.size(), self.size() ) )
		for i, j in product( range( self.size() ), range( self.size() ) ):
			T[ i ][ j ] = ( self.basis[i].T( self.basis[ j ] ) ).sum()
		return T

	def vmatrix( self ):
		V = mtrx( ( self.size(), self.size() ) )
		for i, j in product( range( self.size() ), range( self.size() ) ):
			V[ i ][ j ] = ( self.basis[i].Vp( self.basis[ j ], 50 ) ).sum()
		return V

	def hmatrix( self ):
		H = mtrx( ( self.size(), self.size() ) )
		for i, j in product( range( self.size() ), range( self.size() ) ):
			H[ i ][ j ] = ( self.basis[i].T( self.basis[ j ] ) + \
				        self.basis[i].Vp( self.basis[ j ], 50 ) ).sum()
		return H

	def rhomatrix( self ):
		rho = mtrx( ( self.size(), self.size() ) )
		for i, j in product( range( self.size() ), range( self.size() ) ):
			rho[ i ][ j ] = np.conj( self.basis[i].D ) * \
						 self.basis[j].D
		return rho


	def s1matrix( self ):
		S10 = mtrx( ( self.Ndim * self.size(), self.size() ) )
		for i, j in product( range( 0, self.size() * self.Ndim, self.Ndim ), range( self.size() ) ):
			for k in range( self.Ndim ):
				S10[ i + k ][ j ] = ( self.basis[ i ].der1( self.basis[ j ] ) ).call(k)
		return S10

	def  s11matrix( self ):
		S11 = mtrx( ( self.Ndim * self.size(), self.Ndim * self.size() ) )
		for i, j in product( range( 0, self.size() * self.Ndim, self.Ndim ), range( 0, self.size() * self.Ndim, self.Ndim ) ):
			for k, m in product( range( self.Ndim ), range( self.Ndim ) ):
				S11[ i + k ][ j + m ] = self.basis[ i ].der11( self.basis[ j ] )[k][m]

		return S11

	def h1matrix( self ):
		H10 = mtrx( ( self.Ndim * self.size(), self.size() ) )
		for i, j in product( range( 0, self.size() * self.Ndim, self.Ndim ), range( self.size() ) ):
			for k in range( self.Ndim ):
				H10[ i + k ][ j ] = ( self.basis[ i ].H1( self.basis[ j ] ) ).call(k)
		return H10

	def xmatrix( self ):
		X = self.s11matrix() - np.dot( self.s1matrix(), np.dot( sinvert( self.smatrix() ), np.transpose( self.s1matrix() ) ) )
		rho = self.rhomatrix()
		for i, j in product( range( self.size() ), range( self.size() ) ):
			for k, m in product( range( self.Ndim ), range( self.Ndim ) ):
				X[ i * self.Ndim + k ][ j * self.Ndim + m ] = \
					rho[ i ][ j ] * X[ i * self.Ndim + k ][ j * self.Ndim + m ]
		return X
	
	def ymatrix( self ):
		tmp = self.h1matrix() - np.dot( self.s1matrix(), np.dot( sinvert( self.smatrix() ), self.hmatrix() ) )
		Y = mtrx( self.size() * self.Ndim )
		rho = self.rhomatrix()
		for i, j in product( range( self.size() ), range( self.Ndim ) ):
			for k in range( self.size() ):
				Y[ i * self.Ndim + j ] += rho[ i ][ k ] * tmp[ i * self.Ndim + j ][ k ]
		return Y

	def lambda_dot( self ):
		return -1j * np.dot( sinvert( self.xmatrix() ),  self.ymatrix() )

	def qp_dot( self ):
		ldot = self.lambda_dot()
		qdot = lmap( lambda i: ldot[ i ].real / self.basis[ i // self.Ndim ].omega.call( i % self.Ndim ), \
			     range( self.size() * self.Ndim ) )
		pdot = lmap( lambda i: ldot[ i ].imag, range( self.size() * self.Ndim ) )
		return qdot, pdot

	def taumatrix( self ):
		tau = mtrx( ( self.size(), self.size() ) )
		s1 = self.s1matrix()
		ldot = self.lambda_dot()
		for i, j in product( range( self.size() ), range( self.size() ) ):
			for k in range( self.Ndim ):
				tau[i][j] += s1[ i ][ j * self.Ndim + k ] * ldot[ j * self.Ndim + k ]
		return tau

	def vec_c( self ):
		c = lmap( lambda i: self.basis[i].D, range( self.size() ) )
		return c

	def c_dot( self ):
		return -1j * np.dot( sinvert( self.smatrix() ), np.dot( ( self.hmatrix() - 1j * self.taumatrix() ), self.vec_c() ) )

	def create_tmp_basis( self, y_c, y_q, y_p ):
		b_tmp = []
		for i in range( self.size() ):
			q = []
			p = []
			for j in range( self.Ndim ):
				q.append( y_q[ i * self.Ndim + j ] )
				p.append( y_p[ i * self.Ndim + j ] )
			b_tmp.append( gwp( y_c[i], q, p ) )

		self.basis_tmp = basis( 0.0, b_tmp )
#		print( self.basis_tmp)

	def dynamics( self, delta_t ):
		k1_c_dot = self.c_dot()
		k1_q_dot, k1_p_dot = self.qp_dot()

		y_c = lmap( lambda y, k1: y.D + 0.5 * k1 * delta_t, self.basis, k1_c_dot )
		y_q = [ self.basis[i].q.coords[j] + 0.5 * k1_q_dot[ i * self.Ndim + j ] * delta_t  \
		         for i, j in product( range( self.size() ), range( self.Ndim ) ) ]
		y_p = [ self.basis[i].p.coords[j] + 0.5 * k1_p_dot[ i * self.Ndim + j ] * delta_t  \
			 for i, j in product( range( self.size() ), range( self.Ndim ) ) ]

#		print('1__________________')
		self.create_tmp_basis( y_c, y_q, y_p )
#		print('___________________')

		k2_c_dot = self.basis_tmp.c_dot()
		k2_q_dot, k2_p_dot = self.basis_tmp.qp_dot()

		y_c = lmap( lambda y, k2: y.D + 0.5 * k2 * delta_t, self.basis, k2_c_dot )
		y_q = [ self.basis[i].q.coords[j] + 0.5 * k2_q_dot[ i * self.Ndim + j ] * delta_t  \
			for i, j in product( range( self.size() ), range( self.Ndim ) ) ]
		y_p = [ self.basis[i].p.coords[j] + 0.5 * k2_p_dot[ i * self.Ndim + j ] * delta_t  \
			for i, j in product( range( self.size() ), range( self.Ndim ) ) ]

#		print('2__________________')
		self.create_tmp_basis( y_c, y_q, y_p )
#		print('___________________')

		k3_c_dot = self.basis_tmp.c_dot()
		k3_q_dot, k3_p_dot = self.basis_tmp.qp_dot()

		y_c = lmap( lambda y, k3: y.D + k3 * delta_t, self.basis, k3_c_dot )
		y_q = [ self.basis[i].q.coords[j] + k3_q_dot[ i * self.Ndim + j ] * delta_t  \
			for i, j in product( range( self.size() ), range( self.Ndim ) ) ]
		y_p = [ self.basis[i].p.coords[j] + k3_p_dot[ i * self.Ndim + j ] * delta_t  \
			for i, j in product( range( self.size() ), range( self.Ndim ) ) ]

#		print('3__________________')
		self.create_tmp_basis( y_c, y_q, y_p )
#		print('___________________')

		k4_c_dot = self.basis_tmp.c_dot()
		k4_q_dot, k4_p_dot = self.basis_tmp.qp_dot()

		for i in range( self.size() ):
			self.basis[ i ].D += ( k1_c_dot[i] + 2 * k2_c_dot[i] + 2 * k3_c_dot[i] + k4_c_dot[i] ) * delta_t / 6
			for j in range( self.Ndim ):
				self.basis[ i ].q.coords[ j ] += ( k1_q_dot[ i * self.Ndim + j ] + \
								   2 * k2_q_dot[ i * self.Ndim + j ] + \
								   2 * k3_q_dot[ i * self.Ndim + j ] + \
								   k4_q_dot[ i * self.Ndim + j ] ) * delta_t / 6

				self.basis[ i ].p.coords[ j ] += ( k1_p_dot[ i * self.Ndim + j ] + \
								   2 * k2_p_dot[ i * self.Ndim + j ] + \
								   2 * k3_p_dot[ i * self.Ndim + j ] + \
								   k4_p_dot[ i * self.Ndim + j ] ) * delta_t / 6

	def norm( self ):
		return np.dot( np.conj( self.vec_c() ), np.dot( self.smatrix(), self.vec_c() ) )

	def energy( self ):
		return np.dot( np.conj( self.vec_c() ), np.dot( self.hmatrix(), self.vec_c() ) )

	def get_states( self ):
		x, v = svd( self.smatrix(), self.hmatrix() )
		return x, v

	def create_start( self, i = 0 ):
#		x, v = self.get_states()
		N = np.sqrt( self.norm() )
		for j in range( self.size() ):
#			self.basis[j].D = v[j][i]
			self.basis[j].D /= N 
		

	def __str__ ( self ):
		return reduce( lambda a, x: a + str( x ) + '\n', self.basis, '' )


'''
	Functional style

	def smatrix( self ):
		S = matrix( lambda x, y: ( x * y ).mult() )
		return S( self.basis, self.basis )

	def hmatrix( self ):
		H = matrix( lambda x, y: ( x.T(y) + x.Vf(y, 50) ).sum() )
		return H( self.basis, self.basis )
'''


