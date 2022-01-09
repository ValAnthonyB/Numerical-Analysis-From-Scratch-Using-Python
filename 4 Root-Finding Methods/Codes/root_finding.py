import numpy as np

def bisection_method(f, x0, x1, kmax=100000, tol=1.e-12):
    '''Evaluates the roots of the function f using the bisection method.
    Parameters
        f - function
        
        x0 - lower bound
            float
        
        x1 - upper bound 
            float
        
    Returns
        x2 - root
    
    '''
    f0 = f(x0)
    
    for k in range(0,kmax):
        # Set the midpoint
        x2 = (x0 + x1) / 2
        
        # Evaluate f(x_2)
        f2 = f(x2)
        
        # If f0 * f2 is negative then the root exists at (x0, x2)
        if f0 * f2 < 0:
            x1 = x2
            
        # If f0 * f2 is positive then the root exists at (x2, x1). We search there.
        else:
            x0 = x2
            f0 = f2 # Since we are searching at (x2, x1), the lower bound must be f2
            
            
        # Stopping criteria
        x2new = (x0 + x1) / 2 # check if the next x2new is almost equal to x2. If not do not proceed, the loop repeats
        if np.abs((x2new - x2) / x2new) < tol:
            # If condition is satisfied, then x2 is x*
            break
        
    print("Solution found after", k, "iterations.")
    return x2new


def secant_method(f, x0, x1, kmax=int(1e6), tol=1.e-12):
    '''Provides the approximate root using the secant method.
    
    Parameters
    -----------
        f: function
                Function that we want to find the root.
                
        x0: number
                First initial guess.
                
        x1: number
                Second initial guess.
                
        kmax: integer
                Maximum iterations.
        
        tol: integer
                stopping criteria
        
        
    Returns
    --------
        x: number
            root
    
    '''
    
    for k in range(0, kmax):
        # Secant method
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        
        
        # Stopping criteria
        if np.abs((x2 - x1) / x2) < tol:
            break
            
        else:
            x0 = x1
            x1 = x2
            
    print("Solution found after", k, "iterations.")
    return x2



def riddlers_method(f, x0, x1, kmax=int(1e6), tol=1.e-12):
    '''Provides the approximate root using Riddlers' method.
    
    Parameters
    -----------
        f: function
                Function that we want to find the root.
                
        x0: number
                First endpoint of the bracketed root.
                
        x1: number
                Last endpoint of the bracketed root.
                
        kmax: integer
                Maximum iterations.
        
        tol: integer
                stopping criteria
        
        
    Returns
    --------
        x: number
            root
    
    '''
    
    for k in range(1, kmax):
        # Find the midpoint
        x2 = (x0 + x1) / 2
        
        # Riddlers' method
        x3 = x2 + (x1 - x2)* (f(x2) / f(x0)) / np.sqrt((f(x2) / f(x0))**2 - (f(x1) / f(x0)) )
        
        
        # Stopping criteria
        if np.abs((x3 - x2) / x3) < tol:
            break
            
        else:
            # Change the bracket (x0, x1)_{k+1} = (x2, x3)_{k}
            x0 = x2 
            x1 = x3
            
    print("Solution found after", k, "iterations.")
    return x3