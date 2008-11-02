"""Generate a diffusion stencil

Supports isotropic diffusion (FE,FD), anisotropic diffusion (FE, FD), and
rotated anisotropic diffusion (FD).

The stencils include redundancy to maintain readability for simple cases (e.g.
isotropic diffusion).

-div Q A Q grad u

Q = [cos(theta) -sin(theta)]
    [sin(theta)  cos(theta)]

A = [1          0        ]
    [0          eps      ]

"""

from numpy import array, zeros
from scipy import cos, sin

__docformat__ = "restructuredtext en"

__all__ = ['diffusion_stencil_2d']

def diffusion_stencil_2d(epsilon=1.0, theta=0.0, type='FE'):
    """
    Parameters
    ----------
    type : {'FE','FD'}
        Specifies the discretization as Q1 finite element (FE) or 2nd order
        finite difference (FD)
    dim : {1,2,3}
        The dimension of the problem.  2D or 3D
    epsilon  : float, optional
        Anisotropic diffusion coefficient: -div A grad u, 
        where A = [1 0; 0 `epsilon`].  The default is isotropic, `epsilon`=1.0
    theta : float, optional
        Rotation angle `theta` in radians defines -div Q A Q^T grad,
        where Q = [cos(`theta`) -sin(`theta`); sin(`theta`) cos(`theta`)].  
        The default is `theta` = 0.0

    Return
    ------
    stencil : numpy array
        A 3x3 diffusion stencil

    See Also
    --------
    stencil_grid, poisson

    Notes
    -----
    Not all combinations are supported.

    Examples
    --------
    >>>> from diffusion_stencil import diffusion_stencil
    >>>> sten = diffusion_stencil('FD',dim=2,epsilon=0.0001,theta=pi/6)
    >>>> print sten

    """
    
    eps   = float(epsilon) #for brevity
    theta = float(theta)

    if(type=='FE'):
        """FE approximation to:

            - (eps c^2 +     s^2) u_xx + 
            -2(eps - 1) c s       u_xy + 
            - (    c^2 + eps s^2) u_yy

        [ -c^2*eps-s^2+3*c*s*(eps-1)-c^2-s^2*eps,        2*c^2*eps+2*s^2-4*c^2-4*s^2*eps, -c^2*eps-s^2-3*c*s*(eps-1)-c^2-s^2*eps]
        [       -4*c^2*eps-4*s^2+2*c^2+2*s^2*eps,        8*c^2*eps+8*s^2+8*c^2+8*s^2*eps,       -4*c^2*eps-4*s^2+2*c^2+2*s^2*eps]
        [ -c^2*eps-s^2-3*c*s*(eps-1)-c^2-s^2*eps,        2*c^2*eps+2*s^2-4*c^2-4*s^2*eps, -c^2*eps-s^2+3*c*s*(eps-1)-c^2-s^2*eps]

        c = cos(theta)
        s = sin(theta)
        """
        C = cos(theta)
        S = sin(theta)
        CS = C*S
        CC = C**2
        SS = S**2

        a =  (-1*eps - 1)*CC + (-1*eps - 1)*SS + ( 3*eps - 3)*CS
        b =  ( 2*eps - 4)*CC + (-4*eps + 2)*SS
        c =  (-1*eps - 1)*CC + (-1*eps - 1)*SS + (-3*eps + 3)*CS
        d =  (-4*eps + 2)*CC + ( 2*eps - 4)*SS
        e =  ( 8*eps + 8)*CC + ( 8*eps + 8)*SS
        
        stencil = array([[a,b,c],
                         [d,e,d],
                         [c,b,a]]) / 6.0

    elif type == 'FD':
        """FD approximation to:

        - (eps c^2 +     s^2) u_xx +
        -2(eps - 1) c s       u_xy +
        - (    c^2 + eps s^2) u_yy

          c = cos(theta)
          s = sin(theta)

        A = [ 1/2(eps - 1) c s    -(c^2 + eps s^2)    -1/2(eps - 1) c s  ]
            [                                                            ]
            [ -(eps c^2 + s^2)       2 (eps + 1)    -(eps c^2 + s^2)     ]
            [                                                            ]
            [  -1/2(eps - 1) c s    -(c^2 + eps s^2)  1/2(eps - 1) c s   ]
        """
        C = cos(theta)
        S = sin(theta)
        CC = C**2
        SS = S**2
        CS = C*S
        
        a =  0.5*(eps - 1)*CS
        b = -(eps*SS + CC)
        c = -a
        d = -(eps*CC + SS)
        e =  2.0*(eps + 1)

        stencil = array([[a,b,c],
                         [d,e,d],
                         [c,b,a]])

    return stencil
