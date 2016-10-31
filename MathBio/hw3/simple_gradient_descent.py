# BASED ON CODE FROM: http://en.wikipedia.org/wiki/Gradient_descent

#   f(x)=x**4-3*x**3+2 , with derivative f'(x)=4*x**3-9*x**2

import sys

if len(sys.argv) > 1:
   gamma = float(sys.argv[1])
else:
   # default step size
   gamma = 0.01

# From calculation, we expect that the local minimum occurs at x=9/4

xOld = 0
xNew = 6 # The algorithm starts at x=6

precision = 0.00001
 
def f_prime(x):
    return 4 * x**3 - 9 * x**2
 
istep = 0
while abs(xNew - xOld) > precision:
    xOld = xNew
    grad = f_prime(xNew)
    xNew = xOld - gamma * grad
    istep += 1
    print istep, "xOld,xNew,grad:", xOld,xNew,grad

print("Local minimum occurs at ", xNew)
