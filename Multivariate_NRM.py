'''
 Program to solve system of coupled non-linear equations

 3*x1 - cos(x2*x3) - 3/2 = 0
 4*x1^2 - 625*x2^2 + 2x3 -1 = 0
 20*x3 + exp(-x1*x2) + 9 = 0 

 '''

import numpy as np 
from numpy.linalg import inv

def f1(x1, x2, x3):
	return (3*x1) - np.cos(x2*x3) - 1.5

def f2(x1, x2, x3):
	return (4*pow(x1,2)) - (625*pow(x2,2)) + 2*x3 - 1

def f3(x1, x2, x3):
	return (20*x3) + np.exp(-x1*x2) + 9

def jac(x1, x2, x3):
	# Defining the Jacobian Matrix

	J = np.ones((3,3))

	# Row 1
	J[0,0] = 3
	J[0,1] = -np.sin(x2*x3)*x3
	J[0,2] = -np.sin(x2*x3)*x2

	# Row 2
	J[1,0] = 8*x1
	J[1,1] = -625*2*x2
	J[1,2] = 2

	# Row 3
	J[2,0] = -np.exp(-x1*x2)*x2;
	J[2,1] = -np.exp(-x1*x2)*x1;
	J[2,2] = 20

	return J

# Initial Guess Values

x1 = 1
x2 = 0
x3 = 0

# Initial Error

err = 9e9

# Old Guess Values Column vector

X_old = np.ones((3,1))
X_old[0] = x1
X_old[1] = x2
X_old[2] = x3

# Main Functions Column Vector
F = np.copy(X_old)

alpha = 1.1
tol = 1e-16

# Newton Rhapson loop
iter = 1

while err > tol:
	J = jac(X_old[0], X_old[1], X_old[2])

	#Computing the main functions
	F[0] = f1(X_old[0], X_old[1], X_old[2])
	F[1] = f2(X_old[0], X_old[1], X_old[2])
	F[2] = f3(X_old[0], X_old[1], X_old[2])

	#Computing New Values
	X_new = X_old - alpha*np.matmul(inv(J), F)

	# Computing Max absolute error
	err = np.max(np.abs(X_new - X_old))

	# Updating old values 
	X_old = X_new

	# log message
	log_message = 'Iteration # = {0} x1 = {1} x2 = {2} x3 = {3}'.format(iter, X_new[0], X_new[1], X_new[2])
	print (log_message)

	iter = iter + 1

# Final Results
x1_fin = X_new[0]
x2_fin = X_new[1]
x3_fin = X_new[2]