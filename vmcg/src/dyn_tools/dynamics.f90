SUBROUTINE odeint( NumG, Ndim, Step, c0, cdot, l0, ldot, c, l )

	INTEGER, INTENT(IN) :: NumG, Ndim
	DOUBLE PRECISION, INTENT(IN) :: Step
	DOUBLE COMPLEX, INTENT(IN) :: c0( NumG ), l0( NumG * Ndim ), &
				      cdot( NumG ), ldot( NumG * Ndim )
	DOUBLE COMPLEX, INTENT(OUT) :: c( NumG ), l( NumG * Ndim )

	EXTERNAL :: F
	INTEGER :: ISTATE, ITASK, IPAR, IOPT, &
		   MF, N, LZW, LRW, LIW, LRP
	DOUBLE PRECISION :: T, TOUT, RTOL, ATOL
	INTEGER, ALLOCATABLE :: IWORK(:)
	DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
	DOUBLE COMPLEX, ALLOCATABLE :: ZWORK(:), RPAR(:), Y(:)

	T = 0.0D0
	TOUT = Step
	ITOL = 1
	RTOL = 1.0D-14
	ATOL = 1.0D-14
	ITASK = 1
	ISTATE = 1
	IOPT = 0
	MF = 20
	N =  ( Ndim + 1 ) * NumG 
	LZW = 8 * N
	LRW = 20 + N
	LIW = 30
	LRP = N**2

	ALLOCATE( Y(N), RPAR(:), ZWORK(LZW), RWORK(LRW), IWORK(LIW) )

	ZWORK(:) = (0.0D0, 0.0D0)
	RWORK(:) = 0.0D0
	IWORK(:) = 0

	Y(1:NumG) = c0(:)
	Y(NumG+1:N) = l0(:)
	RPAR(1:NumG) = cdot(:)
	RPAR(NumG+1:N) = ldot(:)

	CALL ZVODE(F, NumG, Y, T, TOUT, ITOL, RTOL, &
		   ATOL, ITASK, ISTATE, IOPT, &
		   ZWORK, LZW, RWORK, LRW, IWORK, LIW, &
		   DUMMY, MF, RPAR, IPAR)

	c(:) = Y(1:NumG)
	l(:) = Y(NumG+1:N)

	RETURN

END SUBROUTINE odeint

SUBROUTINE F( N, T, Y, YDOT, RPAR, IPAR )

	INTEGER, INTENT(IN) :: N, IPAR
	DOUBLE PRECISION, INTENT(IN) :: T
	DOUBLE COMPLEX :: Y(N), YDOT(N), RPAR(N)

	YDOT(:) = RPAR(:)

END SUBROUTINE F
