SUBROUTINE mean_v( NPTS, Ndim, xi1, eta1, omega1, xi2, eta2, omega2, V )

	INTEGER, INTENT(IN) :: NPTS, Ndim
	DOUBLE PRECISION, INTENT(IN) :: omega1(Ndim), omega2(Ndim)
	DOUBLE COMPLEX, INTENT(IN) :: xi1(Ndim), xi2(Ndim), eta1(Ndim), eta2(Ndim)
	DOUBLE COMPLEX, INTENT(OUT) :: V(Ndim)

	INTEGER :: i
	DOUBLE PRECISION :: x(NPTS), x_corr(NPTS,Ndim), wts(NPTS), lb(Ndim), rb(Ndim)
	DOUBLE COMPLEX :: wp1(NPTS,Ndim), wp2(NPTS,Ndim)

	CALL p_quadrature_rule(npts,x,wts)

	V(:) = DCMPLX( 0.0D0, 0.0D0 )
	DO i = 1, Ndim

		lb(i) = MIN( DBLE( xi1(i) ) / omega1(i) - 5.0 / omega1(i), &
			     DBLE( xi2(i) ) / omega2(i) - 5.0 / omega2(i) )
		rb(i) = MAX( DBLE( xi1(i) ) / omega1(i) + 5.0 / omega1(i), &
			     DBLE( xi2(i) ) / omega2(i) + 5.0 / omega2(i) )

		x_corr(:,i) = shift( NPTS, x, lb(i), rb(i) )

		wp1(:,i) = CONJG( wave_packet( npts, x_corr(:,i), xi1(i), eta1(i), omega1(i) ) )
		wp2(:,i) = wave_packet( npts, x_corr(:,i), xi2(i), eta2(i), omega2(i) )

		V(i) = V(i) + SUM( wts(:) * wp1(:,i) * potential( NPTS, x_corr(:,i) ) * wp2(:,i) )
		V(i) = V(i) * 0.5 * ( rb(i) - lb(i) )
	END DO

	RETURN

	CONTAINS
		FUNCTION wave_packet( npts, x, ksi, eta, omega )

			INTEGER :: npts
			DOUBLE PRECISION :: x(npts), omega
			DOUBLE COMPLEX :: wave_packet(npts), ksi, eta

			wave_packet(:) = ZEXP( -0.5 * omega * x(:) * x(:) + ksi * x(:) + eta )

		END FUNCTION wave_packet

		FUNCTION potential( NPTS, x )

			INTEGER :: NPTS
			DOUBLE PRECISION :: x(NPTS), potential(NPTS)
			potential(:) = 0.5 * x(:) * x(:)


		END FUNCTION potential

		FUNCTION shift( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift(N)

			shift(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift


END SUBROUTINE mean_v
