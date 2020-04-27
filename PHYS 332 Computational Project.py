# Jimmy Lilly: PHYS 332 Computational Project - Problem 1 (Due 12/6/19)
#Program to evaluate electric field along z-axis of uniformly charged sphere with hole at top
import numpy as np
import matplotlib.pyplot as plt

# Establish arrays to fill in with z and E_z values
z_array= []
E_z_array_numerical = []
E_z_array_analytical = []
result = 0.0

# Declare radius of sphere (Using 1m for convenience)
R = 1.0

# Measure E_z for many z values (z = 0 to z = 5R)
for i in range(0,101):
    
    # Define z as multiple of sphere's radius R
    z = .01*i*R
    z_array.append(z)
    print("Evaluating field at z = "+str(z)+"R")
    
    # Define function to integrate
    def f(theta):
        numerator = np.sin(theta)*(z-(R*np.cos(theta)))
        denominator = (R**2)+(z**2)-(2*R*z*np.cos(theta))
        return numerator*(denominator**-1.5)
    
    # Define trapezoidal integration function (a: lower bound, b: upper bound, n: # of trapezoids)
    def numerical(a,b,n):
        
        # Step size
        h = (b-a)/float(n)
        
        # First and last terms of trapezoidal sum calculated directly
        s = 0.5*(f(a)+f(b))
        
        # Integrate over number of trapezoids
        for theta in range(0,n):
            
            # Evaluate function at one step further from previous step
            s += f(a+(theta*h))
        return(s*h)
    
    # Define analytical solution to integral and evaluate it at each z
    def analytical(Z):
        
        # Define cos(1 degree) in radians
        value = np.cos(np.pi/180)

        # Express components of first term in analytical solution to integral
        numerator1 = ((value*Z) - R)
        denominator1 = (Z**2)*np.sqrt((Z**2)-(2*value*R*Z)+(R**2))
        
        # Express components of second term in analytical solution to integral
        numerator2 = (Z + R)
        denominator2 = ((Z**2)*np.sqrt((Z**2)+(2*R*Z)+(R**2)))
        
        
        # Calculate result
        result = (numerator1/denominator1) + (numerator2/denominator2)
        return result
    
    # Set a, b, and n (Hole is from 0 to 1 degree, so 'a' = 1 degree)
    A = np.pi/180.0
    B = np.pi
    N = 500
    
    # Add result of each integral to an array and print the result
    E_z_array_numerical.append(numerical(A,B,N))
    print("Numerical result = "+str(numerical(A,B,N)))
    
    # Add result of each calculation to an array and print the result
    E_z_array_analytical.append(analytical(z))
    print("Analytical result = "+str(analytical(z)))

# Make scatter plot of E_z vs. z for both the numerical and analytical solutions
plt.title("Electric Field from Uniformly Charged Sphere with Hole")
plt.xlabel("z (as multiple of sphere radius, R)")
plt.ylabel("E(z) (in units of kQ/2 * m^-2)")
plt.scatter(z_array,E_z_array_numerical,label="Numerical Solution")
plt.scatter(z_array,E_z_array_analytical,label="Analytical Solution")
plt.legend(loc="upper right")
plt.show()