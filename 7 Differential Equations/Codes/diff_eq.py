import numpy as np

def rk4(f, r0, t):
    '''Integrates a system of n first-order differential equations using the RK4 method.
    
    Parameters
    -----------
        f: function
                RHS of the system of differential equations.
                
        r0: array
                Array of n initial values.
                
        t: array
                Values of t.
        
    Returns
    --------
        r: array
            Solution
    
    '''
    # n is the number of equations
    n = len(r0)
    
    # N is the number of grid points
    N = len(t)
    
    # Step size
    h = (t[-1] - t[0]) / N
    
    # Setup r. Each row takes the form (xi, yi), etc
    r = np.zeros((len(t), n))
    r[0] = r0 # Change zeroth row to initial condition
        
    # RK4 algorithm
    for i in range(0, N-1):
        k1 = h * f(r[i],  t[i])
        k2 = h * f(r[i] + 0.5*k1, t[i] + 0.5*h)
        k3 = h * f(r[i] + 0.5*k2, t[i] + 0.5*h)
        k4 = h * f(r[i] + k3, t[i] + h)
        r[i+1] = r[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        
    return r.T