# Multivariate Newton-Rhapson Solver

This program is to implement Multivariate Newton Rhapson Solver for the following system of equations:

$\frac{dy_1}{dt} = -0.04y_1 + 10^4 y_2y_3$

$\frac{dy_1}{dt} = 0.04y_1 - 10^4 y_2y_3 - 3 \times 10^7 y^2_2$

$\frac{dy_3}{dt} = 3 \times 10^7 y^2_2$

When observing the given system of equations, we can say that they are coupled and non-linear. We establish that they are coupled because, y1, y2, y3 can only be solved by solving
the three equations together. We say that it is non-linear because the unknown elements in the equations are getting multiplied by themselves or each other.

Here, the Newton Rhapson method is implemented for each time step set for the given simulation time. The difference is that, in Multivariate, the approximations are calculated in a 
matrix form as there are more than one unkowns. The general principle for the Multivariate Newton Rhapson is given below:

   $X_{n+1} = X_n - \alpha*J^{-1} * F$

where 
- $X_n$ is the old value,
- $X_(n+1)$ is the new value which we are guessing with respect to the old value,
- $\alpha$ is the relaxation factor,
- $J^{-1}$ is the inverse Jacobian matrix and 
- $F$ is a matrix which contains the given system of equations.

In the current challenge, we have calculated the Jacobian matrix numerically as analytical derivatives are not always available. The principle of the numerical Jacobian is expressed below:
                                        
 $\frac{df_i}{dx_j} = \frac{f(x_j+h) - f(x_j)}{h}$
 
For the current case, we are running the simulation for 10 minutes. In the Newton Rhapson solver, the relaxation factor is taken as 1.0. We are also using Implicit Euler or otherwise known
as Backward Differencing method for the outer time integration loop. _Backward Differencing method_ is a finite approximation method which uses the function values at x and x − h where h is
the step size.The general principle for this method is mentioned below: -                                                

$f'(x) = \frac{f(x)- f(x-h)}{h}$

Given the step size (h), note that this formula uses the values at and the point at the previous step. As it moves in the backward direction, it is called the backward difference method.
The given initial conditions for the unknowns are $y_1 = 1, y_2 = 0, y_3 = 0$

Here, we are initializing the ODE's in the form $f'(x) = x_{old \ value} + dt*f(x) - x_{new \ value}$
The general layout of the code is given below: -
-	Outer Time Integration Loop
-	Assembling the ODE System
-	Newton Rhapson solver till convergence
-	Get new values and update them
-	End of Outer Loop

