import numpy as np
import matplotlib.pyplot as plt

g = 9.82
l = 1

def Jacobian(alpha, alphadot, dt):
    return lambda alpha, alphadot, dt: np.array([[1, -dt/2],[dt/2/l*g*np.cos(alpha),1]])

def Jacobi_method(J, b, tolerance):
    U = np.triu(J)-np.identity(2)
    L = np.tril(J)-np.identity(2)
    D = np.array([[1,0],[0,1]]) 

    guess = np.array([1,1])
#    oldGuess = guess
    for _ in range(100):
        guess = np.matmul(np.matmul(-D,(L+U)),guess)+ np.matmul(D,b)        
#        if np.linalg.norm(guess - oldGuess, 2) < tolerance*1e-3:
#            break
#        oldGuess = guess
    #print("D:" + str(D))
    #print("L:" + str(L))
    #print("U:" + str(U))
    #print("b:" + str(b))
    #print("guess:" +str(guess))#, np.linalg.solve(J, np.transpose(b)))
    return guess

def Newton_iter(J, initial, F, tolerance):
    delta = Jacobi_method(J, -F, tolerance)
    #delta = np.linalg.solve(J, -F)#to test Newtons method
    return initial + delta
    
def Newtons(J, start, F, tolerance, dt):
    
    old = start

    while True:
        new = Newton_iter(J(old[0],old[1], dt), old, F(old[0],old[1]), tolerance)
        old = new
#        print(new)
        if np.linalg.norm(F(new[0],new[1]),2) < tolerance:
            break
    return new

#Jacobi_method(Jacobian(1,1,0.01))

def pendulum(deltaT, alphaOld, alphaDotOld):
    def F(alpha, alphaDot):
        return np.array([alpha - alphaOld - deltaT/2*(alphaDotOld + alphaDot), 
                         alphaDot - alphaDotOld + deltaT/2*(9.82*(np.sin(alphaOld) + np.sin(alpha)))])#9,82 = g/l
    return F

initialPos = np.array([0, np.pi*2])
pos = initialPos
steps = 200
stop = 10 #seconds
stepsize = stop/steps
tol = 1e-2
posList = [[],[]]
for t in range(steps):
    F = pendulum(stepsize, pos[0], pos[1])
    #FPrim = Jacobian(pos[0], pos[1], stepsize)
    posList[0].append(pos[0])
    posList[1].append(pos[1])
    pos = Newtons(Jacobian(pos[0], pos[1], stepsize), pos, F, tol, stepsize)
#    print(pos)

#print(posList)
plt.figure(dpi=200)
plt.plot(posList[0],posList[1])
plt.xlabel("alpha (rad)")
plt.ylabel("angular velocity (rad/s)")


plt.figure(dpi = 200)
plt.plot(np.linspace(0, stop, steps), posList[0], label = "alpha (rad)")
plt.plot(np.linspace(0, stop, steps), posList[1], label = "alphaDot (rad/s)")
plt.xlabel("t (s)")
plt.legend()