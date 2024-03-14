import numpy as np
import matplotlib.pyplot as plt
import cmath
def matrix(UR,gamma,omega,epsilonk):
    pass

import numpy as np

def GreenFunction(UR, gamma, epsilon, omega):
    # Define the matrix elements
    a11 = (2 * (UR - 2j * gamma + epsilon + omega) * (gamma * (-4j * UR + 5 * gamma) - 2 * (UR + 2j * gamma) * epsilon - epsilon**2 + omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a12 = (2 * (UR - 1j * gamma) * (-5 * gamma**2 + 4j * gamma * epsilon + epsilon**2 + 2 * UR * (2j * gamma + epsilon) - omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a13 = -(8j * gamma * (UR**2 + gamma**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a14 = (8 * gamma * (1j * UR + gamma) * (UR + 2j * gamma + epsilon - omega)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    
    a21 = (2 * (UR + 1j * gamma) * (11 * gamma**2 + 4j * gamma * epsilon + epsilon**2 + 2 * UR * (2j * gamma + epsilon) - omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a22 = (2 * (UR - 2j * gamma + epsilon - omega) * (gamma * (-4j * UR + 5 * gamma) - 2 * (UR + 2j * gamma) * epsilon - epsilon**2 + omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a23 = (8j * (UR + 1j * gamma) * gamma * (UR - 2j * gamma + epsilon - omega)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a24 = -(8j * gamma * (4 * gamma**2 + (UR + epsilon - omega)**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    
    a31 = -(8j * gamma * (4 * gamma**2 + (UR + epsilon + omega)**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a32 = (8 * gamma * (1j * UR + gamma) * (UR + 2j * gamma + epsilon + omega)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a33 = (2 * (UR + 2j * gamma + epsilon + omega) * (-5 * gamma**2 - 4j * gamma * epsilon + epsilon**2 + 2 * UR * (-2j * gamma + epsilon) - omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a34 = (2 * (UR - 1j * gamma) * (4j * UR * gamma - 11 * gamma**2 - 2 * UR * epsilon + 4j * gamma * epsilon - epsilon**2 + omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    
    a41 = (8j * (UR + 1j * gamma) * gamma * (UR - 2j * gamma + epsilon + omega)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a42 = -(8j * gamma * (UR**2 + gamma**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a43 = (2 * (UR + 1j * gamma) * (4j * UR * gamma + 5 * gamma**2 - 2 * UR * epsilon + 4j * gamma * epsilon - epsilon**2 + omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    a44 = (2 * (UR + 2j * gamma + epsilon - omega) * (-5 * gamma**2 - 4j * gamma * epsilon + epsilon**2 + 2 * UR * (-2j * gamma + epsilon) - omega**2)) / ((3 * gamma**2 + epsilon * (2 * UR + epsilon))**2 + 2 * (5 * gamma**2 - epsilon * (2 * UR + epsilon)) * omega**2 + omega**4)
    
    # Construct the matrix
    matrix = np.array([[a11, a12, a13, a14],
                       [a21, a22, a23, a24],
                        [a31, a32, a33, a34],
                        [a41, a42, a43, a44]])

    return matrix


def eigenspectrum(epsilon,UR,gamma):
    rootinside=gamma**2-epsilon*(epsilon+2*UR)
    root = cmath.sqrt(rootinside)
    lastterm = 4*gamma*root
    firstterm=epsilon*(epsilon+2*UR)-5*gamma**2
    return cmath.sqrt(firstterm+lastterm), cmath.sqrt(firstterm-lastterm)


def plotting_imaginary_spectrum(UR,gamma):
    imgs1=[]
    imgs2=[]
    epsilons=[]
    for epsilon in np.linspace(0,3.5,1000):
        a,b=eigenspectrum(epsilon,UR,gamma)
        imgs1.append(a.imag)
        imgs2.append(b.imag)
        epsilons.append(epsilon)
    plt.plot(epsilons,imgs1)
    plt.plot(epsilons,imgs2)
    plt.plot(epsilons,-np.array(imgs1))
    plt.plot(epsilons,-np.array(imgs2))
    plt.legend(["1","2","3","4"])
    plt.show()


def diagonal_Greenfunction_with_gamma(epsilon,omega,UR):
    diags=[]
    gammas=[]
    for gamma in np.linspace(0.5,1.5,100):
        diags.append(np.diag(GreenFunction(epsilon,omega,UR,gamma)))
        gammas.append(gamma)
    plt.plot(gammas,diags)
    plt.legend(["1","2","3","4"])
    plt.show()
