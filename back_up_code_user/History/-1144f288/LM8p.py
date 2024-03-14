import numpy as np
import matplotlib.pyplot as plt

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
