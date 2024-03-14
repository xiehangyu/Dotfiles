def epsilonL(U,a):
    dim=np.shape(a)[0]
    temp_matrix=np.kron(np.eye(dim),a)
    temp_matrix=U.conj().transpose().dot(temp_matrix.dot(U))
    return_matrix=np.zeros((dim,dim),dtype=complex)
    zero_vector=np.zeros(dim,dtype=complex)
    for i in range(dim):
