import numpy as np
import matplotlib.pyplot as plt
from interpolate import interpolate_vec
from collision import has_solution

g = 9.82
l = 1

#Defines the jacobian for our problem
def Jacobian(alpha, alphadot, dt):
    return lambda alpha, alphadot, dt: np.array([[1, -dt/2],[dt/2/l*g*np.cos(alpha),1]])

#Jacobi's method for solving linear equation systems
def Jacobi_method(J, b, tolerance):
    U = np.triu(J)-np.identity(2)
    L = np.tril(J)-np.identity(2)
    D = np.array([[1,0],[0,1]]) 

    #Arbitrary guess, will converge under suitable conditions
    guess = np.array([1,1])

    for _ in range(100):
        guess = np.matmul(np.matmul(-D,(L+U)),guess)+ np.matmul(D,b)        
    return guess

#One iteration of Newton's method for finding roots
def Newton_iter(J, initial, F, tolerance):
    delta = Jacobi_method(J, -F, tolerance)
    return initial + delta
   
#Newton's method for finding roots of F
#Jacobian = J
#Initial value = start
def Newtons(J, start, F, tolerance, dt):
    
    old = start

    while True:
        new = Newton_iter(J(old[0],old[1], dt), old, F(old[0],old[1]), tolerance)
        old = new
        if np.linalg.norm(F(new[0],new[1]),2) < tolerance:
            break
    return new


def pendulum(deltaT, alphaOld, alphaDotOld):
    def F(alpha, alphaDot):
        return np.array([alpha - alphaOld - deltaT/2*(alphaDotOld + alphaDot), 
                         alphaDot - alphaDotOld + deltaT/2*(9.82*(np.sin(alphaOld) + np.sin(alpha)))])
        #9,82 = g/l
    return F



#Intial conditions
initialPos = np.array([0, 2*np.pi])
pos = initialPos
stop = 5 #seconds to stop at
stepsize = 1/500

#Tolerance for Newton's method in every step
tol = 1e-2

posList = [[],[]]

#The pendulum is allowed to pass through the obstacle the frame after hitting it
#As to not get stuck inside it
free_pass = False

#Time of every step
times = [0]

#Simulation loop
while times[-1] < stop:
    times.append(times[-1]+stepsize)
    
    #Update the position
    F = pendulum(stepsize, pos[0], pos[1])
    posList[0].append(pos[0])
    posList[1].append(pos[1])
    pos = Newtons(Jacobian(pos[0], pos[1], stepsize), pos, F, tol, stepsize)

    #Check for collisions during the last update
    if len(posList[0]) > 2:
        #Take the three last points
        current_pos = np.array([[posList[0][-1]],[posList[1][-1]]])
        last_pos = np.array([[posList[0][-2]],[posList[1][-2]]])
        lastlast_pos = np.array([[posList[0][-3]],[posList[1][-3]]])

        #Interpolate them as a function of t
        interp = interpolate_vec(current_pos, last_pos, lastlast_pos, 0, -1, -2)

        alpha_obst = 2*np.pi*np.floor(pos[0]/(2*np.pi)+1/2)-np.pi/6

        #See if it collides with the obstacle within the last timestep
        if has_solution(interp, 0, 1, alpha_obst)[0] and not free_pass:
            #Elastically change the veloicty
            free_pass = True
            times[-1] = times[-2] + has_solution(interp,0,1,alpha_obst)[1]*stepsize
            pos[0] = interp(has_solution(interp,0,1,alpha_obst)[1])[0][0] 
            pos[1] = -interp(has_solution(interp,0,1,alpha_obst)[1])[1][0]
        elif free_pass:
            free_pass = False

#Plotting
plt.figure(dpi=200)
plt.plot(posList[0],posList[1])
plt.xlabel("alpha (rad)")
plt.ylabel("angular velocity (rad/s)")

plt.figure(dpi = 200)
plt.plot(times[1:], posList[0], label = "alpha (rad)")
plt.plot(times[1:], posList[1], label = "alphaDot (rad/s)")
plt.xlabel("t (s)")
plt.legend()
plt.show()
