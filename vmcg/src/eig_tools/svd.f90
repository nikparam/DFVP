SUBROUTINE svd(NumG, S, H, eigen_states, eigen_vectors)

	INTEGER, INTENT(IN) :: NumG
	DOUBLE COMPLEX, INTENT(IN) :: S(NumG,NumG), H(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	DOUBLE COMPLEX :: SI(NumG,NumG), H_prime(NumG,NumG)

	INTEGER :: ok
	DOUBLE COMPLEX :: WORK( 2 * NumG ), DUMMY
	DOUBLE PRECISION :: RWORK( 2 * NumG )

	CALL sinvert(NumG, S, SI)

	H_prime = MATMUL( SI, H )

	CALL ZGEEV('N', 'V', NumG, H_prime, NumG, eigen_states, DUMMY, 1, eigen_vectors, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	RETURN

END SUBROUTINE svd
