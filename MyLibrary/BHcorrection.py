# =============================================================================
# MIT License
# 
# Copyright (c) 2024 Gael Rigaud  <https://www.f08.uni-stuttgart.de/organisation/team/Rigaud/>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# =============================================================================


import numpy as np
from scipy.special import lambertw



def BHcorrection(g_bh,T_E,E,intervalg,absorption_coef,diffusion_coef,tau=0.5,tol=10**(-12)):
    """
    Correction of the Beam-Hardening effect

    Parameters
    ----------
    g_bh : the CT-data affected by beam-hardening (numpy.ndarray())
    T_E, E : the spectrum and the energy range (numpy.ndarray())
    intervalg : a tuple (a,b) which locates the effective area in which
            the approximation is computed.
    absorption_coef : the global absorption coefficient, a constant.
    diffusion_coef : the global diffusion coefficient, a constant.
    tau : parameter to control the approximation on the monochromatic part. 
        tau = 0 means that the Compton cross-section is taken at the lowest 
        energy, while tau = 1 means that the cross-section is taken at the
        highest energy. tau>1 is possible but does make physically much sense.
        The default is 0.5.
    tol : Tolerance for the convergence of the Levenberg-Marquardt algorithm.
        The default is 10**(-12).

    Returns
    -------
    a numpy.ndarray()

    """
    
    dE  = E[1]-E[0]
    iE  = 1/E**3
    C_E   = DiffusionCS(E)
    M = 100
    x = np.linspace(intervalg[0],intervalg[1],M)
    P = np.zeros(M)
    for i in range(M):     
        P[i] = dE*np.sum(T_E*np.exp(-x[i]*iE))

    params = LM(P,x,np.array([10**(-5),1]),tol)    
    alpha = (C_E[0] + tau*(C_E[-1] - C_E[0]))*diffusion_coef
    beta  = params[0]*absorption_coef
    c     = params[1]
    
    print('Parameters of the correction:'+' alpha=',alpha,' beta=', beta, ' c=',c)

    return EvalCorrection(g_bh, alpha, beta, c) 

def DiffusionCS(E):
    """
    The diffusion cross-section.

    Parameters
    ----------
    E : energy (numpy.ndarray)

    Returns
    -------
    numpy.ndarray

    """
    return  ((1+E/511)/(E/511)**2 *(2*(1+E/511)/(1+2*E/511)-np.log(1+2*E/511)/(E/511))+np.log(1+2*E/511)/(2*E/511)-(1+3*E/511)/(1+2*E/511)**2)

def EvalCorrection(g,alpha,beta,c):
    """
    The correction function evaluated on data g.

    Parameters
    ----------
    g : data (numpy.ndarray)
    alpha, beta, c : parameters of the approximation 

    Returns
    -------
    corrected data (numpy.ndarray)

    """
    coef = c*beta/alpha
    return np.real(coef*lambertw(np.exp(1/coef+g/c)/coef) - 1)/beta

def ApproxPoly(x,p):
    """
    Polynom involved in the approximation.
    
    """
    p[p<0] = 0
    return np.power(1+p[0]*x,-p[1])

def grad_ApproxPoly(x,p):
    """
    Gradient of the polynom involved in the approximation.
    
    """
    J = np.zeros((x.shape[0],p.shape[0]))
    b,c = (p[0],p[1])
    J[:,0] = -c*x*(b*x+1)**(-c - 1)
    J[:,1] = -(b*x + 1)**(-c) * np.log(b*x + 1)
    return J


def LM(y,x,p,sigma=100,tol=10**(-10)):    
    """
    The Levenberg-Marquardt algorithm for curve fitting.

    Parameters
    ----------
    y : the signal to approximate (numpy.ndarray)
    x : the coordinates (numpy.ndarray)
    p : the approximation parameters (numpy.ndarray)
    sigma : standard deviation in the chi2 function.
        The default is 100.
    tol : tolerance to convergence.
        The default is 10**(-10).

    Returns
    -------
    p : the approximation parameters at convergence. (numpy.ndarray)

    """
    n = p.shape[0]
    W = 1/sigma**2 * np.eye(y.shape[0])
    d = 2*tol
    i = 0
    mu = 0.1
    b0 = 0.3
    b1 = 0.9

    while np.linalg.norm(d) > tol:
        
        y2 = ApproxPoly(x,p)  
        J = grad_ApproxPoly(x,p)
        Jt= J.T
        JtWJ = Jt @ W @ J
        JtW  = Jt @ W
        
        while True:
            d = np.linalg.solve(JtWJ+mu*np.eye(n), JtW @ (y-y2))
            rho = (chi2(x,p,y,W) - chi2(x,p+d,y,W))/abs(d.dot(mu*d + JtW @ (y-y2)))
            break
            if rho <= b0: mu = 2*mu
            else: break
        p = p + d
        
        if rho > b1: mu = 0.5*mu
        i += 1
        if i>1000: raise Exception('The curve fitting of the polychromatic part is taking too long. Please consider different parameters.')
        
        
    return p

def chi2(x,p,y,W):
    """
    chi2 function associated to the polynom involved in the approximation.
    """
    y2 = ApproxPoly(x,p)
    return (y-y2).dot( W @ (y-y2))
