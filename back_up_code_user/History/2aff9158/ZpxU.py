def epsilonL(U,a):
    dim=np.shape(a)[0]
    temp_matrix=np.kron(np.eye(dim),a)
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    