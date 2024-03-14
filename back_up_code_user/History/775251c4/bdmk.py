import numpy as np
import scipy
from scipy.linalg import expm

def tilde_U(U):
    '''
    temp_U=np.zeros(np.shape(U),dtype=np.complex128)
    temp_U2=U.transpose()
    temp_U[0,0]=temp_U2[0,0]
    temp_U[0,1]=temp_U2[0,1]
    temp_U[0,2]=temp_U2[1,0]
    temp_U[0,3]=temp_U2[1,1]
    temp_U[1,0]=temp_U2[0,2]
    temp_U[1,1]=temp_U2[0,3]
    temp_U[1,2]=temp_U2[1,2]
    temp_U[1,3]=temp_U2[1,3]
    temp_U[2,0]=temp_U2[2,0]
    temp_U[2,1]=temp_U2[2,1]
    temp_U[2,2]=temp_U2[3,0]
    temp_U[2,3]=temp_U2[3,1]
    temp_U[3,0]=temp_U2[2,2]
    temp_U[3,1]=temp_U2[2,3]
    temp_U[3,2]=temp_U2[3,2]
    temp_U[3,3]=temp_U2[3,3]
    return temp_U.transpose()
    '''
    temp_U=np.copy(U)
    temp_U[0,1]=U[2,0]
    temp_U[0,3]=U[2,2]
    temp_U[1,1]=U[3,0]
    temp_U[1,3]=U[3,2]
    temp_U[2,0]=U[0,1]
    temp_U[2,2]=U[0,3]
    temp_U[3,0]=U[1,1]
    temp_U[3,2]=U[1,3]
    return temp_U
def generaltildeU(U):
    dim=int(np.sqrt(np.shape(U)[0]))
    new_matrix=np.zeros(np.shape(U),dtype=np.complex128)
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(dim):
                    new_matrix[j*dim+l,i*dim+k]=U[k*dim+l,i*dim+j]
    return new_matrix

def generalwhether_right_invariant(U):
    dim=int(np.sqrt(np.shape(U)[0]))
    U_tilde=generaltildeU(U)
    right_vector=U_tilde.conjugate().transpose().dot(U_tilde)
    reference_vector=np.copy(right_vector)
    right_vector=np.kron(right_vector,np.eye(dim))
    right_vector=np.dot(np.kron(np.eye(dim),U_tilde.conjugate().transpose()),right_vector)
    right_vector=np.dot(right_vector,np.kron(np.eye(dim),U_tilde))
    return reference_vector, right_vector-np.kron(np.eye(dim),reference_vector)



def whether_right_invariant(U):
    U_tilde=tilde_U(U)
    right_vector=U_tilde.conjugate().transpose().dot(U_tilde)
    reference_vector=np.copy(right_vector)
    right_vector=np.kron(right_vector,np.eye(2))
    right_vector=np.dot(np.kron(np.eye(2),U_tilde.conjugate().transpose()),right_vector)
    right_vector=np.dot(right_vector,np.kron(np.eye(2),U_tilde))
    return reference_vector, right_vector-np.kron(np.eye(2),reference_vector)

def whether_right_invariant_3(U):
    U_tilde=tilde_U(U)
    right_vector=U_tilde.conjugate().transpose().dot(U_tilde)
    right_vector=np.kron(right_vector,np.eye(2))
    right_vector=np.dot(np.kron(np.eye(2),U_tilde.conjugate().transpose()),right_vector)
    right_vector=np.dot(right_vector,np.kron(np.eye(2),U_tilde))
    reference_vector=np.copy(right_vector)
    right_vector=np.kron(right_vector,np.eye(2))
    right_vector=np.dot(np.kron(np.eye(4),U_tilde.conjugate().transpose()),right_vector)
    right_vector=np.dot(right_vector,np.kron(np.eye(4),U_tilde))
    return reference_vector, right_vector-np.kron(np.eye(2),reference_vector)

def generalwhether_left_invariant(U):
    dim=int(np.sqrt(np.shape(U)[0]))
    U_tilde=generaltildeU(U)
    left_vector=U_tilde.dot(U_tilde.conjugate().transpose())
    reference_vector=np.copy(left_vector)
    left_vector=np.kron(left_vector,np.eye(2))
    left_vector=np.dot(np.kron(np.eye(2),U_tilde),left_vector)
    left_vector=np.dot(left_vector,np.kron(np.eye(2),U_tilde.conjugate().transpose()))
    return reference_vector, left_vector-np.kron(np.eye(2),reference_vector)


def whether_left_invariant(U):
    U_tilde=tilde_U(U)
    left_vector=U_tilde.dot(U_tilde.conjugate().transpose())
    reference_vector=np.copy(left_vector)
    left_vector=np.kron(left_vector,np.eye(2))
    left_vector=np.dot(np.kron(np.eye(2),U_tilde),left_vector)
    left_vector=np.dot(left_vector,np.kron(np.eye(2),U_tilde.conjugate().transpose()))
    return reference_vector, left_vector-np.kron(np.eye(2),reference_vector)

def whether_left_invariant_3(U):
    U_tilde=tilde_U(U)
    left_vector=U_tilde.dot(U_tilde.conjugate().transpose())
    left_vector=np.kron(left_vector,np.eye(2))
    left_vector=np.dot(np.kron(np.eye(2),U_tilde),left_vector)
    left_vector=np.dot(left_vector,np.kron(np.eye(2),U_tilde.conjugate().transpose()))
    reference_vector=np.copy(left_vector)
    left_vector=np.kron(left_vector,np.eye(2))
    left_vector=np.dot(np.kron(np.eye(4),U_tilde),left_vector)
    left_vector=np.dot(left_vector,np.kron(np.eye(4),U_tilde.conjugate().transpose()))
    return reference_vector, left_vector-np.kron(np.eye(2),reference_vector)


def overall_testing(U):
    CNOT=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    output1=np.kron(U,np.eye(2)).dot(CNOT)
    output2=np.kron(np.eye(2),U).dot(CNOT)
    input1=CNOT.dot(np.kron(U,np.eye(2)))
    input2=CNOT.dot(np.kron(np.eye(2),U))
    a1,b1=whether_right_invariant(output1)
    print("right invariant output1 {}\n".format(np.trace(np.dot(b1,b1.conjugate().transpose()))))
    a2,b2=whether_right_invariant(output2)
    print("right invariant output2 {}\n".format(np.trace(np.dot(b2,b2.conjugate().transpose()))))
    a3,b3=whether_right_invariant(input1) 
    print("right invariant input1 {}\n".format(np.trace(np.dot(b3,b3.conjugate().transpose()))))
    a4,b4=whether_right_invariant(input2)
    print("right invariant input2 {}\n".format(np.trace(np.dot(b4,b4.conjugate().transpose()))))
    a5,b5=whether_left_invariant(output1)
    print("left invariant output1 {}\n".format(np.trace(np.dot(b5,b5.conjugate().transpose()))))
    a6,b6=whether_left_invariant(output2)
    print("left invariant output2 {}\n".format(np.trace(np.dot(b6,b6.conjugate().transpose()))))
    a7,b7=whether_left_invariant(input1)
    print("left invariant input1 {}\n".format(np.trace(np.dot(b7,b7.conjugate().transpose()))))
    a8,b8=whether_left_invariant(input2)
    print("left invariant input2 {}\n".format(np.trace(np.dot(b8,b8.conjugate().transpose()))))

def translation_operator(n):
    '''
    1,2,3,n->n,1,2,3
    '''
    matrix_size=(1<<n)
    translation_operator_matrix=np.zeros((matrix_size,matrix_size))
    for i in range(matrix_size):
        after_translation=(i>>1)+(i&1)*(1<<(n-1))
        translation_operator_matrix[after_translation,i]=1
    return translation_operator_matrix

def classification_of_ab():
    X=np.array([[0,1],[1,0]])
    Z=np.array([[1,0],[0,-1]])
    CNOT=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
    fp1=open("./plot_data/classification_ab_second_left.txt","w")
    fp2=open("./plot_data/classification_ab_second_right.txt","w")
    gp1=open("./plot_data/classification_ab_third_left.txt","w")
    gp2=open("./plot_data/classification_ab_third_right.txt","w")
    for a in np.linspace(0,2*np.pi,500):
        for b in np.linspace(0,2*np.pi,500):
            print("a={},b={}".format(a,b))
            useful_matrix=np.kron(np.array([[0,0],[0,1]]),expm(-1j*(np.cos(b)*X+np.sin(b)*Z)*a))+np.kron(np.array([[1,0],[0,0]]),np.eye(2))
       #     temp_matrix=np.kron(np.eye(2),expm(-1j*(np.cos(b)*X+np.sin(b)*Z)*a))
       #     useful_matrix=temp_matrix.conj().transpose().dot(CNOT.dot(temp_matrix))
            temp1,temp2=whether_left_invariant(useful_matrix)
            if (np.trace(np.dot(temp2,temp2.conj().transpose()))/np.trace(np.dot(temp1,temp1.conj().transpose()))<1e-7):
                fp1.write("{}\t{}\n".format(a,b))
            temp1,temp2=whether_right_invariant(useful_matrix)
            if (np.trace(np.dot(temp2,temp2.conj().transpose()))/np.trace(np.dot(temp1,temp1.conj().transpose()))<1e-7):
                fp2.write("{}\t{}\n".format(a,b))
            temp1,temp2=whether_left_invariant_3(useful_matrix)
            if (np.trace(np.dot(temp2,temp2.conj().transpose()))/np.trace(np.dot(temp1,temp1.conj().transpose()))<1e-7):
                gp1.write("{}\t{}\n".format(a,b))
            temp1,temp2=whether_right_invariant_3(useful_matrix)
            if (np.trace(np.dot(temp2,temp2.conj().transpose()))/np.trace(np.dot(temp1,temp1.conj().transpose()))<1e-7):
                gp2.write("{}\t{}\n".format(a,b))
    fp1.close()
    fp2.close()
    gp1.close()
    gp2.close()
            

def eigenvalues_of_concatenated_CNOT(n):
    CNOT=np.array([[1.0,0,0,0],[0,1.0,0,0],[0,0,0,1.0],[0,0,1.0,0]])
    CNOT=tilde_U(CNOT)
    translation_operator_matrix=translation_operator(2*n)
    CNOT_matrix=np.copy(CNOT)
    for i in range(n-1):
        CNOT_matrix=np.kron(CNOT_matrix,CNOT)
    CNOT_matrix2=translation_operator_matrix.dot(CNOT_matrix.dot(translation_operator_matrix.transpose()))
    CNOT_return_matrix=np.dot(CNOT_matrix2,CNOT_matrix)
    return np.linalg.eigvals(np.dot(CNOT_return_matrix.transpose(),CNOT_return_matrix))


def printklrv(D,alpha,beta,gamma):
    '''
    Calculate S1, where S1 is the element which makes the first brackets non zero.
    '''
    print("(k,l)=")
    for k in range(D):
        for l in range(D):
            if(((2*alpha*k+beta*l)%D==0) and ((beta*k+2*gamma*l)%D==0)):
                print("({},{})\t".format(k,l),end='')
    print("\n\n")
    print("(r,v)=")
    '''
    Calculate S2, where S2 is the element which makes the second brackets non zero.
    '''
    my_set=set()
    for k in range(D):
        for l in range(D):
            r=beta*k+2*gamma*l
            v=-2*alpha*k-beta*l
            r=r%D
            v=v%D
            my_set.add((r,v))
    my_str='\t'.join(str(x) for x in my_set)
    print(my_str)
    print("\n\n")
    print("(m,n)=")
    '''
    Calculate S3, where S3 is the element which makes the off diagonal elements non zero.
    '''
    for k in range(D): 
        for l in range(D): 
            if(((2*gamma*l+beta*k+k)%D==0) and ((2*alpha*k+beta*l-l)%D==0)):
                print("({},{})\t".format(k,l),end='')

def search_non_vanishing_offdiagonal_element(D):
    alpha=0
    beta=0
    gamma=0
    if D%2==0:
        step=1/2
    else:
        step=1
    while alpha<=D:
        beta=0
        while beta<=D:
            gamma=0
            while gamma<=D:
                key2=0
                for k in range(D):
                    for l in range(D):
                        if(((2*alpha*k+beta*l)%D==0) and ((beta*k+2*gamma*l)%D==0)):
                            if(k!=0 or l!=0):
                                key2=1
                                break
                    if key2==1:
                        break
                if key2==0:
                    gamma+=step
                    continue
                key=0
                for k in range(D): 
                    for l in range(D):
                        if(((2*gamma*l+beta*k+k)%D==0) and ((2*alpha*k+beta*l-l)%D==0)):  
                            if k!=0 or l!=0:
                                key=1
                                break
                    if key==1:
                        break
                if key==1:
                    print("alpha={},beta={},gamma={}".format(alpha,beta,gamma))
                gamma+=step
            beta+=1
        alpha+=step
    
                             
def klrv(D,alpha,beta,gamma):
    '''
    Calculate S1, where S1 is the element which makes the first brackets non zero.
    '''
    my_set1=set()
    my_set2=set()
    for k in range(D):
        for l in range(D):
            if(((2*alpha*k+beta*l)%D==0) and ((beta*k+2*gamma*l)%D==0)):
                if k!=0 or l!=0:
                    my_set1.add((k,l))
    '''
    Calculate S2, where S2 is the element which makes the second brackets non zero.
    '''
    my_set3=set()
    for k in range(D):
        for l in range(D):
            r=beta*k+2*gamma*l
            v=-2*alpha*k-beta*l
            r=r%D
            v=v%D
            if r!=0 or v!=0:
                my_set2.add((r,v))
    '''
    Calculate S3, where S3 is the element which makes the off diagonal elements non zero.
    '''
    for k in range(D): 
        for l in range(D): 
            if(((2*gamma*l+beta*k+k)%D==0) and ((2*alpha*k+beta*l-l)%D==0)):
                if k!=0 or l!=0:
                    my_set3.add((k,l))
    return my_set1, my_set2, my_set3

def condition1(D):
    if D%2==0:
        step=1/2
    else:
        step=1
    alpha=0
    while alpha<D:
        beta=0
        while beta<D:
            gamma=0
            while gamma<D:
                my_set1, my_set2, my_set3=klrv(D,alpha,beta,gamma)
                if len(my_set1)==0:
                    gamma+=step
                    continue
                if len(my_set3)==0:
                    gamma+=step
                    continue
                key=0
                for temp in my_set1:
                    if temp not in my_set2:
                        key=1
                if key==0:
                    gamma+=step
                    continue
                print("alpha={},beta={},gamma={}".format(alpha,beta,gamma))
                gamma+=step
            beta+=1
        alpha+=step
    