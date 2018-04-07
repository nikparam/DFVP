SUBROUTINE cholesky(NumG, S, H, eigen_states, eigen_vectors)

	INTEGER, INTENT(IN) :: NumG
	DOUBLE COMPLEX, INTENT(IN) :: S(NumG,NumG), H(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	DOUBLE COMPLEX :: tmpS(NumG, NumG), H_prime(NumG,NumG), U_prime(NumG,NumG), DUMMY

	INTEGER :: ok, INFO
	DOUBLE PRECISION :: RWORK(2*NumG)
	DOUBLE COMPLEX :: WORK(2*NumG)

	tmpS(:,:) = S(:,:)
	CALL ZPOTRF( 'L', NumG, tmpS, NumG, INFO )
	DO i = 1, NumG
		DO j = i+1, NumG
			tmpS(i,j) = 0.0
		END DO
	END DO

	INFO = 0
	CALL ZTRTRI( 'L', 'N', NumG, tmpS, NumG, INFO )

	DO i = 1, NumG
		DO j = i+1, NumG
			tmpS(i,j) = 0.0
		END DO
	END DO

	H_prime = MATMUL( tmpS , MATMUL( H, TRANSPOSE( CONJG( tmpS ) ) ) )

	CALL ZGEEV('N', 'V', NumG, H_prime, NumG, eigen_states, DUMMY, 1, U_prime, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	eigen_vectors = MATMUL( TRANSPOSE( CONJG( tmpS ) ), U_prime )

	RETURN

END SUBROUTINE cholesky
