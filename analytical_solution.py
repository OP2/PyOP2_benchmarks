from numpy import cos, exp, matrix, pi, sqrt

def advection_diffusion(x, t):
    dx= (matrix(x)-matrix((0.49 + t,0.5)))
    r=sqrt(dx*dx.T)
    D = 0.1 # Diffusivity
    A = 0.1 # Normalisation
    return A*(exp((-r**2)/(4*D*t))/(4*pi*D*t))
