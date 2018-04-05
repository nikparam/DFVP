SUBROUTINE mean_v( NPTS, q1, p1, omega1, q2, p2, omega2, V )

	INTEGER, INTENT(IN) :: NPTS
	DOUBLE PRECISION, INTENT(IN) :: q1(3), q2(3), p1(3), p2(3), omega1, omega2
	DOUBLE PRECISION, INTENT(OUT) :: V

	INTEGER :: i
	DOUBLE PRECISION :: x(NPTS), x_corr(NPTS,3), wts(NPTS), lb(3), rb(3)

	CALL p_quadrature_rule(npts,x,wts)

	i = 0
20 	CONTINUE 	
		i = i + 1
		lb(i) = MIN( q1(i) - 5.0/omega1, q2(i) - 5.0/omega2 )
		rb(i) = MAX( q1(i) + 5.0/omega1, q2(i) + 5.0/omega2 )

		x_corr(:,i) = shift( NPTS, x, lb(i), rb(i) )
	IF ( i < 3 ) GOTO 20

	V = 0.

	DO i = 1, NPTS
		DO j = 1, NPTS
			DO k = 1, NPTS
				V = V + 0.5 * ()
			END DO
		END DO
	END DO


	RETURN

	CONTAINS
		FUNCTION wave_packet( npts, x, ksi, eta, omega )

			INTEGER :: npts
			DOUBLE PRECISION :: x(npts), omega
			DOUBLE COMPLEX :: wave_packet(npts), ksi, eta

			wave_packet(:) = ZEXP( -0.5 * omega * x(:) * x(:) + ksi * x(:) + eta )

		END FUNCTION wave_packet

		FUNCTION potential_map(npts, x, params)

			INTEGER :: npts, i
			DOUBLE PRECISION :: x(npts), params(15)
			DOUBLE PRECISION :: potential_map(npts)

			DO i = 1, npts			
				CALL potential_energy( x(i), params, potential_map(i) )
			END DO

		END FUNCTION potential_map

		FUNCTION shift( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift(N)

			shift(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift


END SUBROUTINE mean_v
