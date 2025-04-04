import math
import numpy as np

#def dscb(float(x), float(alphaL), float(alphaR), float(nL), float(nR), float(N), float(mean), float(sigma)):
def dscb(x,*par):
    #print(x[0])
    #print(par)
    #t = (x[0] - par[5])/par[6]
    result = []
    t = (x - par[5])/par[6]
    #print(t)
    A = par[2]/par[0]
    B = par[3]/par[1]
    N = par[4]
    for i in t:
        if(-par[0] <= i <= par[1]):
            f = math.exp(-0.5*i**2)
        elif(i < -par[0]):
            f = math.exp(-0.5*par[0]**2)*(1/A*(A - par[0] -i))**(-par[2])
        elif(i > par[1]):
            f = math.exp(-0.5*par[1]**2)*(1/B*(B - par[1] +i))**(-par[3])
        else:
            f = np.inf
        result.append(N*f)
    #result = N*f
    if (type(result) is complex):
        result = (result).real
    #print(result)
    return result

def splusb(x, bcurve, bscale, *par):
    fbkg = bcurve.interpolate(x) * bscale
    result = []
    t = (x - par[5])/par[6]
    A = par[2]/par[0]
    B = par[3]/par[1]
    N = par[4]
    for n,i in enumerate(t):
        if(-par[0] <= i <= par[1]):
            f = math.exp(-0.5*i**2)
        elif(i < -par[0]):
            f = math.exp(-0.5*par[0]**2)*(1/A*(A - par[0] -i))**(-par[2])
        elif(i > par[1]):
            f = math.exp(-0.5*par[1]**2)*(1/B*(B - par[1] +i))**(-par[3])
        else:
            f = np.inf
        result.append(N*f + fbkg[n])
    #result = N*f
    #if (type(result) is complex):
    #    result = (result).real
    #print(result)
    return result

"""
def dscb(x,alphaL,alphaH,nL,nH,N,mean,sigma):
    t = (x - mean)/sigma
    A = nL/alphaL
    B = nH/alphaH
    if(-alphaL <= t <= alphaH):
        f = N*exp(-0.5*t**2)
    elif(t < -alphaL):
        f = N*exp(alphaL**2)(1/A*(A - alphaL -t))**(-nL)
    elif(t > alphaH):
        f = N*exp(alphaH**2)(1/B*(B - alphaH +t))**(-nH)
        
    return f
"""

def crystalball(x, *par):
    """
    Single-sided Crystal Ball function
    Parameters:
        x: float or array-like, input value(s)
        par[0]: a (tail transition parameter)
        par[1]: n (power-law exponent for the tail)
        par[2]: N (normalization factor)
        par[3]: mean (Gaussian center)
        par[4]: sigma (Gaussian width)
    """
    result = []
    t = (x - par[3]) / par[4]
    a = par[0]
    n = par[1]
    N = par[2]
    
    for i in t:
        if i > -a:
            f = math.exp(-0.5 * i**2)  # Gaussian core
        else:
            A = (n / abs(a))
            B = (n / abs(a)) - abs(a)
            f = math.exp(-0.5 * a**2) * (1 / (A * (B - i)))**n  # Power-law tail
        
        result.append(N * f)
    
    #return np.array(result) if isinstance(x, np.ndarray) else result[0]
    return result
'''
def gausTanh(x, bcurve, bscale, *par):
    norm = par[0]
    mean = par[1]
    sigma = par[2]
    gaus_ = norm*np.exp(np.exp(-((x - mu) ** 2) / (2 * sigma ** 2)))
    tanh_ = np.tanh(x)
    fbkg = bcurve.interpolate(x) * bscale
    return (gaus_ * tanh_) + fbkg
'''
 
def parabola(p, *pars):  # a, pmin, cmin
    #print(*pars)
    return pars[0]*(p - pars[1])**2 + pars[2]

def sqrt_func(x,a):
    return a*np.sqrt(x)