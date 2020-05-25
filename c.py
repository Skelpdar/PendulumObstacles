import numpy as np

g = 9.82
l = 1

def Jacobian(alpha, alphadot, dt):
    return np.array([[1, dt/2],[dt/2/l*g*np.cos(alpha),1]])

def Jacobi_method(J, b):
    U = np.triu(J)-np.identity(2)
    L = np.tril(J)-np.identity(2)
    D = np.array([[1,0],[0,1]]) 

    guess = np.array([[0.01],[0.01]])

    while error > tolerance:
        guess = np.matmul(np.matmul(-D,(L+U)),guess)+ np.matmul(D,b)
        error = ?

    return guess

def Newton_iter(J, initial, F):
    delta = Jacobi_method(J, -F)
    return initial + delta
    
def Newtons(J, start, F, tolerance, dt):
    
    old = start

    while True:
        new = Newton_iter(Jacobian(old[0][0],old[1][0], dt), old, F(old[0][0],old[1][0]))

        if np.linalg.norm(new - old, 2) < tolerance:
            break

Jacobi_method(Jacobian(1,1,0.01))
