#Newton Rhapson Method
# x^3 + 4*x^2 +2X - 5

import matplotlib.pyplot as plt 

def f(x):
	return pow(x,3) + 4*pow(x,2) + 2*x - 5

def f1(x):
	return 3*pow(x,2) + 8*x + 2

# Newtin Rhapson Iteration

x_guess = 80
alpha = 1.5
tol = 1e-6
iter = 1

while (abs(f(x_guess))>tol):
	x_guess = x_guess - alpha*(f(x_guess)/f1(x_guess))
	iter = iter + 1

print(iter)
print(x_guess)