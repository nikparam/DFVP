SUBROUTINE lowdin(NumG, S, H, eigen_states, eigen_vectors)

	DOUBLE COMPLEX, INTENT(IN) :: S(NumG,NumG), H(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	DOUBLE COMPLEX :: S_mhalf(NumG,NumG), s_eigen(NumG), TR(NumG, NumG), &
			  H_prime(NumG,NumG), H_dprime(NumG,NumG), &
			  U_prime(NumG,NumG), DUMMY

	INTEGER :: ok
	DOUBLE PRECISION :: RWORK(2*NumG)
	DOUBLE COMPLEX :: WORK(2*NumG)

	ok = 0
	CALL ZGEEV('N', 'V', NumG, S, NumG, s_eigen, DUMMY, 1, TR, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	S_mhalf(1:NumG,1:NumG) = (0.0D0, 0.0D0)
	DO i = 1, NumG
		S_mhalf(i,i) = 1.0D0 / DSQRT( DBLE( s_eigen(i) ) )
	END DO

	H_prime = MATMUL( TRANSPOSE( CONJG( TR ) ), MATMUL( H, TR ) )
	H_dprime = MATMUL( S_mhalf, MATMUL( H_prime, S_mhalf ) )

	WORK = ( 0.0D0, 0.0D0 )
	RWORK = 0.0D0
	DUMMY = ( 0.0D0, 0.0D0 )
	ok = 0

	CALL ZGEEV('N', 'V', NumG, H_dprime, NumG, eigen_states, DUMMY, 1, U_prime, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	eigen_vectors = MATMUL( TR, MATMUL( S_mhalf, U_prime ) )

	RETURN 

END SUBROUTINE lowdin
