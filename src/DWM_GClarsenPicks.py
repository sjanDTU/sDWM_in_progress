"""GC Larsen wake model applied to offshore wind farms (WindFarm object)
@moduleauthor:: Juan P. Murcia <jumu@dtu.dk>
References:
[1] Larsen GC. "A simple stationary semi-analytical wake model", 2009
"""

import numpy as np
def get_Rw_NOJ(x,D,CT):
    kj=0.05
    RW=np.zeros((len(x),1))
    Area=pi*D*D/4.0
    a=(1.0-(np.sqrt(1.0-CT)))/2.0
    k=np.sqrt((1.0-a)/(1.0-2.0*a))
    for i in np.arange(0,len(x),1):
        if (x[i]<= 0):
            RW[i] = 0.0
        else:
          RW[i]=0.5*k*D + kj*x[i]

    return RW, kj

def gaussN(R, func, varargin, NG=4):
    """Calculate numerically the gauss integration.
    [1] eq. 38

    Inputs
    ----------
    R (float): Wind turbine radius [m]
    func (function): Wind speed function
    varargin: Other arguments for the function besides [r,te]
    NG (int): Number of Ga

    Outputs
    ----------
    Ua (float):
    """
    A = np.pi*R**2
    #coefficients
    if  NG==4: #for speed give the full values
        rt = np.array([[ -0.339981043584856, -0.861136311594053,
            0.339981043584856, 0.861136311594053]])
        te = rt.T
        w = np.array([[0.652145154862546, 0.347854845137454,
            0.652145154862546, 0.347854845137454]])
    else:
        rt,w = np.polynomial.legendre.leggauss(NG)
        rt = np.array([rt])
        te = rt.T
        w = np.array([w])

    return (np.pi/4.0)*(R**2./A)*w*w.T*func(R*(rt+1.0)/2.0,
        np.pi*(te+1.0),*varargin)*(rt+1.0)


    return R96

def get_Rw(x,R,TI,CT,pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake radius at a location

    Inputs
    ----------
    x (float): Distance between turbines and wake location in the wind direction
    R (float): Wind turbine radius
    TI (float): Ambient turbulence intensity
    CT (float): Outputs WindTurbine object's thrust coefficient

    Outputs
    ----------
    Rw (float): Wake radius at a location
    """
    _ones = np.ones(np.shape(x))
    D = 2*R                          # D is the diameter of the wake generating rotor
    Area=np.pi*D*D/4.0               # Area is the rotor area

    #CT=4.0*a*(1.-a)

    m=1.0/(np.sqrt(1.0-CT))
    k=np.sqrt((m+1.0)/2.0)

    R96 = get_R96(R,CT,TI,pars)      # The wake radius at 9.6D down stream distance is (empirally) approximated with the following expression

    x0=(9.6*D)/((2.0*R96/(k*D))**3.0-1.0)
    term1=(k*D/2.0)**2.5
    term2=(105.0/(2.0*np.pi))**-0.5
    term3=(CT*Area*x0)**(-5.0/6.0)
    c1=term1*term2*term3

    Rw=((105.0*c1*c1/(2.0*np.pi))**0.2)*(CT*Area*(x+x0*_ones))**(1.0/3.0)
    #print 'Rw: '
    #print Rw

    if type(x)==float and x+x0 <= 0.:
        Rw = 0
    elif type(x)==np.ndarray:
        Rw[x+x0*_ones <= 0.] = 0.
    #print 'Rw (before return): '
    #print Rw
    return Rw

def get_R96(R,CT,TI,pars=[0.435449861,0.797853685,-0.124807893,0.136821858,15.6298,1.0]):
    """Computes the wake radius at 9.6D downstream location of a turbine

    Inputs
    ----------
    R (float): Wind turbine radius
    CT (float): Outputs WindTurbine object's thrust coefficient
    TI (float): Ambient turbulence intensity

    Outputs
    ----------
    R96 (float): Wake radius at 9.6D downstream location
    """
    D = 2.0*R

    a1=pars[0]
    a2=pars[1]
    a3=pars[2]
    a4=pars[3]
    b1=pars[4]
    b2=pars[5]
    R96=a1*(np.exp(a2*CT*CT+a3*CT+a4*CT**0.))*(b1*TI+b2)*D

    return R96

def Ua(r,te,zc,us,z0):
    """Function of undisturbed inflow wind speed - log law.

    Inputs
    ----------
    r  (np.array float): Radial coord [m]
    te (np.array float): Angle coord where 0 is the horizontal right [rad]
    zc (float): Height to hub center [m]
    us (float): Friction velocity [m/s]
    z0 (float): Roughness height [m]

    Outputs
    ----------
    Ua  (np.array float): Axial wind speed [m/s]
    """
    kappa = 0.4 # Kappa: von karman constant
    return us / kappa * np.log((zc + r * np.sin(te)) / z0)

def Ua_shear(r,te,zc,uH,alpha):
    """Function of undisturbed inflow wind speed - power law.

    Inputs
    ----------
    r  (np.array float): Radial coord [m]
    te (np.array float): Angle coord where 0 is the horizontal right [rad]
    zc (float): Height to hub center [m]
    uH (float): Wind Speed at Hub Height [m/s]
    alpha (float): Shear Coefficient [-]

    Outputs
    ----------
    Ua  (np.array float): Axial wind speed [m/s]
    """

    return uH * ((zc + r * np.sin(te)) / zc)**alpha
