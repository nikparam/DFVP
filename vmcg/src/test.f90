PROGRAM test

	DOUBLE PRECISION :: V, q1(3), p1(3), q2(3), p2(3)
	DOUBLE COMPLEX :: xi1(3), eta1(3), xi2(3), eta2(3)

	q1 = (/ 0.0D0, 0.0D0, 0.0D0 /)
	p1 = (/ 0.0D0, 0.0D0, 0.0D0 /)
	q2 = (/ 0.0D0, 0.0D0, 0.0D0 /)
	p2 = (/ 0.0D0, 0.0D0, 0.0D0 /)

	xi1(:) = DCMPLX( q1(:), p1(:) )
	xi2(:) = DCMPLX( q2(:), p2(:) )

	phase = 0.25 * DLOG( 0.25 / DATAN( 1.0D0 ) )
	eta1(:) = DCMPLX( phase - 0.5 * q1(:) * q1(:), -q1(:) * p1(:) )
	eta2(:) = DCMPLX( phase - 0.5 * q2(:) * q2(:), -q2(:) * p2(:) )

	CALL mean_v( 48, xi1, eta1, 1.0D0, xi2 , eta2, 1.0D0, V)

	WRITE(*,*) V

END PROGRAM test
