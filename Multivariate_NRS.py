'''
Program to solve system of coupled non-linear equations using Multivariate Newton Rhapson Method
-0.04*y1 + 10^4*y2*y3 = 0
0.04*y1 - 10^4*y2*y3 - 3*10^7*y2^2 = 0
3*10^7*y2^2 = 0

'''

import numpy as np 
from numpy.linalg import inv
import matplotlib.pyplot as plt 

#Step size initialization
dt = 0.6
t_start = 0
t_end = 600  # 10 minutes = 60*10

# Initial Guess Values
y1 = 1
y2 = 0
y3 = 0

def f1(y1_old, y2_old, y3_old, y1, y2, y3, dt):
	return y1_old + dt*(-0.04*y1 + pow(10,4)*y2*y3) - y1

def f2(y1_old, y2_old, y3_old, y1, y2, y3, dt):
	return y2_old + dt*(0.04*y1 - pow(10,4)*y2*y3 - 3*pow(10,7)*pow(y2,2)) - y2

def f3(y1_old, y2_old, y3_old, y1, y2, y3, dt):
	return y3_old + dt*(3*pow(10,7)*pow(y2,2)) - y3

def jac(y1_old,y2_old,y3_old,y1, y2, y3, dt):
	# Defining the Jacobian Matrix

	J = np.ones((3,3))
	d = 1e-8

	# Row 1
	J[0,0] = (f1(y1_old, y2_old, y3_old, y1+d, y2, y3, dt) - f1(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d
	J[0,1] = (f1(y1_old, y2_old, y3_old, y1, y2+d, y3, dt) - f1(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d
	J[0,2] = (f1(y1_old, y2_old, y3_old, y1, y2, y3+d, dt) - f1(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d

	# Row 2
	J[1,0] = (f2(y1_old, y2_old, y3_old, y1+d, y2, y3, dt) - f2(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d
	J[1,1] = (f2(y1_old, y2_old, y3_old, y1, y2+d, y3, dt) - f2(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d
	J[1,2] = (f2(y1_old, y2_old, y3_old, y1, y2, y3+d, dt) - f2(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d

	# Row 3
	J[2,0] = (f3(y1_old, y2_old, y3_old, y1+d, y2, y3, dt) - f3(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d
	J[2,1] = (f3(y1_old, y2_old, y3_old, y1, y2+d, y3, dt) - f3(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d
	J[2,2] = (f3(y1_old, y2_old, y3_old, y1, y2, y3+d, dt) - f3(y1_old, y2_old, y3_old, y1, y2, y3, dt))/d

	return J


def Newton_Rhapson(y1, y2, y3, y1_guess, y2_guess, y3_guess, dt):

	# Old Guess Values Column vector
	Y_old = np.ones((3,1))
	Y_old[0] = y1_guess
	Y_old[1] = y2_guess
	Y_old[2] = y3_guess

	F = np.ones((3,1))

	#Parameters for Newton Raphson solver
	error = 9e9  #Initializing error
	alpha = 1
	tol = 1e-9
	
	iter = 1

	while (error > tol):

		J = jac(y1, y2, y3, Y_old[0],Y_old[1],Y_old[2], dt)

		#Computing the main functions
		F[0] = f1(y1, y2, y3, Y_old[0],Y_old[1],Y_old[2], dt)
		F[1] = f2(y1, y2, y3, Y_old[0],Y_old[1],Y_old[2], dt)
		F[2] = f3(y1, y2, y3, Y_old[0],Y_old[1],Y_old[2], dt)

		# Computing New Values
		Y_new = Y_old - alpha*(np.matmul(inv(J), F))

		# Computing Max absolute error
		error = np.max(np.abs(Y_new - Y_old))

		# Updating old values 
		Y_old = Y_new

		# log message
		log_message = 'Iteration # = {0} y1 = {1} y2 = {2} y3 = {3}'.format(iter, Y_new[0], Y_new[1], Y_new[2])
		print (log_message)

		iter = iter + 1

	return [Y_new[0], Y_new[1], Y_new[2]]


def implicit_euler(y1, y2, y3, t_start, t_end, dt):
	
	# Newton Rhapson Outer Integration Loop
	t = np.arange(t_start,t_end,dt)
	time = len(t)

	Y_1 = np.zeros(time)
	Y_2 = np.zeros(time)
	Y_3 = np.zeros(time)
	
	Y_1[0] = y1
	Y_2[0] = y2
	Y_3[0] = y3

	y1_guess = 0
	y2_guess = 0
	y3_guess = 0

	for i in range(1,time):

		Y_1[i], Y_2[i], Y_3[i] = Newton_Rhapson(Y_1[i-1], Y_2[i-1], Y_3[i-1],y1_guess,y2_guess,y3_guess,dt)

		# Final Results
		y1_guess = Y_1[i]
		y2_guess = Y_2[i]
		y3_guess = Y_3[i]

	return [t, Y_1, Y_2, Y_3]

t, y1, y2, y3 = implicit_euler(1,0,0,0,600,0.6)

plt.plot(t,y1)
plt.plot(t,y2)
plt.plot(t,y3)
plt.xlabel('Simulation Time')
plt.legend(['y1 graph','y2 graph','y3 graph'])
plt.show()


#Y[i+1] = Y[i] + h*f'(Y[i]) Explicit
#Y[i] = Y[i-1] + h*f'(Y[i]) --> -Y[i] + Y[i-1] + h*f'(Y[i]) = 0 Implicit