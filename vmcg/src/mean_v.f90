SUBROUTINE mean_v( NPTS, xi1, eta1, omega1, xi2, eta2, omega2, V )

	INTEGER, INTENT(IN) :: NPTS
	DOUBLE PRECISION, INTENT(IN) :: omega1, omega2
	DOUBLE COMPLEX, INTENT(IN) :: xi1(3), xi2(3), eta1(3), eta2(3)
	DOUBLE COMPLEX, INTENT(OUT) :: V

	INTEGER :: i
	DOUBLE PRECISION :: x(NPTS), x_corr(NPTS,3), wts(NPTS), lb(3), rb(3)
	DOUBLE COMPLEX :: wp1(NPTS,3), wp2(NPTS,3)

	CALL p_quadrature_rule(npts,x,wts)

	i = 0
20 	CONTINUE 	
		i = i + 1
		lb(i) = MIN( DBLE( xi1(i) )/omega1 - 5.0/omega1, DBLE( xi2(i) )/omega2 - 5.0/omega2 )
		rb(i) = MAX( DBLE( xi1(i) )/omega1 + 5.0/omega1, DBLE( xi2(i) )/omega2 + 5.0/omega2 )

		x_corr(:,i) = shift( NPTS, x, lb(i), rb(i) )
	IF ( i < 3 ) GOTO 20

	DO i = 1, 3
		wp1(:,i) = wave_packet(npts, x_corr(:,i), xi1(i), eta1(i), omega1 )
		wp2(:,i) = CONJG( wave_packet(npts, x_corr(:,i), xi2(i), eta2(i), omega2 ) )
	END DO

	V = 0.
	DO i = 1, NPTS
		DO j = 1, NPTS
			DO k = 1, NPTS
				V = V + wts(i) * wts(j) * wts(k) * &
					wp2(i,1) * wp2(j,2) * wp2(k,3) * &
					potential( (/ x_corr(i,1), x_corr(j,2), x_corr(k,3) /) ) * &
					wp1(i,1) * wp1(j,2) * wp1(k,3)
			END DO
		END DO
	END DO

	V = V * 0.125 * ( rb(1) - lb(1) ) * &
			( rb(2) - lb(2) ) * &
			( rb(3) - lb(3) )

	RETURN

	CONTAINS
		FUNCTION wave_packet( npts, x, ksi, eta, omega )

			INTEGER :: npts
			DOUBLE PRECISION :: x(npts), omega
			DOUBLE COMPLEX :: wave_packet(npts), ksi, eta

			wave_packet(:) = ZEXP( -0.5 * omega * x(:) * x(:) + &
							      ksi * x(:) + eta )

		END FUNCTION wave_packet

		FUNCTION potential( x )

			DOUBLE PRECISION :: x( 3 ), potential
			potential = 0.5 * SUM( x(:) * x(:) )


		END FUNCTION potential

		FUNCTION shift( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift(N)

			shift(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift


END SUBROUTINE mean_v
