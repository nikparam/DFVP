SUBROUTINE sinvert( NumG, S, SI )

	INTEGER, INTENT(IN) :: NumG
	DOUBLE cOMPLEX, INTENT(IN) :: S(NumG, NumG )
	DOUBLE COMPLEX, INTENT(OUT) :: SI(NumG, NumG)
	DOUBLE PRECISIOn :: Sd(NumG)
	DOUBLE COMPLEX :: Sdm(NumG,NumG), &
			  U(NumG,NumG), VT(NumG,NumG)

	INTEGER :: i, LW, LIW, LRW, astatus
	INTEGER, ALLOCATABLE :: IWORK(:)
	DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
	DOUBLE COMPLEX, ALLOCATABLE :: WORK(:)

	LW = NumG*( NumG + 3 )
	LIW = 8 * NumG
	LRW = NumG * ( 5 * NumG + 7 )
	ALLOCATE(WORK(LW), IWORK(LIW), RWORK(LRW), STAT=astatus)
	CALL ZGESDD( 'A', NumG, NumG, S, NumG, Sd, U, NumG, VT, NumG, WORK, LW, RWORK, IWORK, INFO )

	Sdm(1:NumG,1:NumG) = (0.0D0, 0.0D0)
	DO i = 1, NumG
		IF ( ABS( Sd(i) ) .LT. ( 1.0D-8 ) ) THEN
			Sdm(i,i) = Sd(i) + 1.0D-8 * EXP( -Sd(i) * 1.0D8)
		END IF
			Sdm(i,i) = 1.0D0 / Sd(i)
	END DO

	SI = MATMUL( TRANSPOSE( CONJG( VT ) ), MATMUL( Sdm, TRANSPOSE( CONJG( U ) ) ) )


END suBROUTINE sinvert
