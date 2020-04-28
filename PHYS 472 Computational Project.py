"""
Author: Jimmy Lilly (www.github.com/jlilly364)

Program Objective: Numerical Integration of the Time-Independent Schrodinger 
                   Equation for a 1D Linear Half-Potential
"""
####
# V(x)=alpha*x for x>0, infinite elsewhere
# E_3 should be 8.6217 eV
# Function is even so: wavefunction initially 0, derivate initally 1.0
####

# Import relevant modules
import matplotlib.pyplot as plt
import numpy as np

# Define functional form of time-independent potential: V(x)
def Potential(x,alpha,exponent):
    # Inputs:
    #    x = value to evaluate potential at
    #    alpha = multiplicative constant for potential
    #    exponent = potential's exponential dependence on x
    # Returns:
    #    value of potential at given x value
    function = alpha*(x**exponent)
    return function

def Test(x):
    # Function to test 'Solver' where V(x)=0
    return 0*x

# Function to solve TISE
def Solver(stepsize,upper,energy,even):
    # Inputs:
    #     stepsize = size of deltaX for iterating loop (in angstroms)
    #     upper = upper limit for integration (in angstroms)
    #     energy = trial energy
    #     even = boolean for if wavefunction is even (True) or odd (False)
    # Returns:
    #     xrange = x values at which Psi was evaluated
    #     Psi_array = wavefunction evaluated for all values in xrange
    
    # Establish initial conditions by testing if function is odd or even
    if even == True: # Even
        Psi = 1.0
        dPsi = 0.0    
        print('This wavefunction is even')
    else: # Odd
        Psi = 0.0
        dPsi = 1.0
        print('This wavefunction is odd')
    
    # 2m/hbar^2 in natural units
    beta = 0.26246 # in (ev*angstroms^2)^-1
    
    # Set up array to add Psi values too
    Psi_array = []
    
    # Establish list of x values
    xrange = np.arange(0,upper+stepsize,stepsize)
    
    # Calculate how many loops will run
    numLoops = int(upper/stepsize)
    
    # Loop to calculate important quantities
    for i in range(0,len(xrange)):
        
        # Tell user which loop is running
        print('Running loop {0} of {1}'.format(i,numLoops))
        
        # Calc. 2nd derivative of Psi using
        d2Psi = beta*(Potential(xrange[i],1.0,1) - energy)*Psi
        
        # What to add to wavefunction
        extra = ((dPsi*stepsize) + ((d2Psi*(stepsize**2))/2))
        
        # Calc. new value of dPsi: dPsi(x0+deltaX)
        dPsi += d2Psi*stepsize
        
        # Calc. new value of Psi: Psi(x0+deltaX) & add to array
        Psi += extra
        Psi_array.append(Psi)
    
    # Add initial value of Psi to array of Psi values
    np.insert(Psi_array,1,Psi)
    
    return xrange,Psi_array

# Function to find general locations of eigenenergies
def BigFishing(central,buffer=1.00,step=0.25): 
    
    # Make list of energies to run through
    energies = np.arange(central-buffer,central+buffer+step,step)
    x = [[] for i in range(len(energies))]
    y = [[] for i in range(len(energies))]
    
    # Set values to plot for each test energy
    for i in range(0,len(energies)):
        x[i],y[i] = Solver(.001,15.0,energies[i],True)

    # Plot wavefunctions for each test energy
    for i in range(len(x)):
        plt.plot(x[i],y[i],label='Energy={0} eV'.format(energies[i]))
        plt.xlim(0,15)
        plt.ylim(-1.5,1.5)
        plt.hlines(0,0,15)
        plt.legend()
        plt.xlabel(r'Distance ($\AA$)')
        plt.ylabel(r'$\psi$(x)')
        plt.title('Wavefunction for Linear Potential')
    
    plt.show()
    
# Function to find more precise eigenenergies
def MediumFishing(central,buffer=.25,step=0.125):
    
    # Make list of energies to run through
    energies = np.arange(central-buffer,central+buffer+step,step)
    x = [[] for i in range(len(energies))]
    y = [[] for i in range(len(energies))]
    
    # Set values to plot for each test energy
    for i in range(0,len(energies)):
        x[i],y[i] = Solver(.001,15.0,energies[i],True)

    # Plot wavefunctions for each test energy
    for i in range(len(x)):
        plt.plot(x[i],y[i],label='Energy={0} eV'.format(np.around(energies[i],4)))
        plt.xlim(0,15)
        plt.ylim(-1.5,1.5)
        plt.hlines(0,0,15)
        plt.xlabel(r'Distance ($\AA$)')
        plt.ylabel(r'$\psi$(x)')
        plt.legend()
    plt.show()
    
# Function to find more precise eigenenergies
def SmallFishing(central,buffer=.001,step=0.0005):
    
    # Make list of energies to run through
    energies = np.arange(central-buffer,central+buffer+step,step)
    x = [[] for i in range(len(energies))]
    y = [[] for i in range(len(energies))]
    
    # Set values to plot for each test energy
    for i in range(0,len(energies)):
        x[i],y[i] = Solver(.001,15.0,energies[i],True)

    # Plot wavefunctions for each test energy
    for i in range(len(x)):
        plt.plot(x[i],y[i],label='Energy={0} eV'.format(np.around(energies[i],4)))
        plt.xlim(0,15)
        plt.ylim(-1,1)
        plt.hlines(0,0,15)
        plt.xlabel(r'Distance ($\AA$)')
        plt.ylabel(r'$\psi$(x)')
        plt.legend()
    plt.show()

#BigFishing(8.0)
#MediumFishing(7.56)
#SmallFishing(7.5275)

# Plot of final eigenenergy wavefunctions
def FinalPlot(save=False):
    E1 = 1.5907
    E2 = 3.6515
    E3 = 7.5275
    E4 = 8.6217
    x1,y1 = Solver(0.001,15.0,1.5907,True)
    x2,y2 = Solver(0.001,15.0,3.6515,False)
    x3,y3 = Solver(0.001,15.0,E4,True)
    x4,y4 = Solver(0.001,15.0,E4,False)
    plt.plot(x1,y1,label=r'E$_1$ = {0} eV'.format(E1))
    plt.plot(x2,y2,label=r'E$_2$ = {0} eV'.format(E2))
    plt.plot(x3,y3,label=r'E$_3$ (even) = {0} eV'.format(E4))
    plt.plot(x4,y4,label=r'E$_3$ (odd) = {0} eV'.format(E4))
    plt.legend()
    plt.xlim(0,15)
    plt.ylim(-1.5,1.5)
    plt.hlines(0,0,15,linestyle='dashed')
    plt.xlabel(r'Distance ($\AA$)')
    plt.ylabel(r'$\psi$(x)')
    plt.title('Wavefunction for Linear Potential')
    if save == True:
        plt.savefig('C:/Users/Jimmy/Physics-Programs/Linear Half Potential Energies')
    plt.show()

FinalPlot(True)