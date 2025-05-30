'''
                                     .   .                 .:
                                     ^^:.:               .^!?~
                                     !^:^:             .^~:.:.  .
                                  :^:.^::!:           .:~^.^:::^!.
                                  .:^::^::!~.:~~^::::.:!~.~Y77::^
                                   ~::^^~!!JJY5YYY55P5?J7!J?~:...
                                   ^^:^^!?J?JJJ?JJYPPPYYY5?~~^~~!~
                                    .::^~??YY?777J???YJ?JYY?7~^::
                                   .:!7!^^JY??!77??J777?YJJ5^..
                                .^!!7??7~::?JJ7!777?777?JJ?J^
                           .:~!!JYYYJYYY?7!!7JYYYYY5YJJ!~^^~!~~~!7?7:
          ::        ...::~!7??!?YJJ555YJJ77????JYY?~.        .:?J!^:^
          .~!7!7777777????77!!~!!~!???JYYPPPPP55!.             ^ .^
            .^!????777777777!~~~~~~!!!!~~!~~~~^~~~:.
               .^!??77????7777?7!~!~^.           :^!77^
                  ...:~!77!^:^7?!^                  ^~?!
                              !?~:                    .:
                             :77

AxolotE Python version. Here, the recurrence relationships are implemented.
'''
from sympy import *
from sympy.series.sequences import RecursiveSeq
import functools

@functools.cache
def ECoeff(l=0, m=0, lp=0, mp=0, t=0, u=0, v=0,
    p=symbols("p"), Dx=symbols("Dx"), Dy=symbols("Dy"), Dz=symbols("Dz"),
    Dr2=symbols("Dr2"), Dxp=symbols("Dxp"), Dyp=symbols("Dyp"),
    Dzp=symbols("Dzp"), Drp2=symbols("Drp2")
):
    """Generates the expressions for the expansion coefficients for the
    product of real solid harmonics (rsh) as a lineal combination of
    Hermite polynomials (Hp), in the form they are found on gaussian
    type functions (gtf).

    Solid Harmonics in the Shift phase convention, which basically
    omits the Condon-Shortley phase coefficient, is used. Unnormalized
    functions are considered.

    If no parameters are passed the coefficient for the product of
    two s real solid harmonics (l=0, m=0, lp=0 and mp=0) is returned.

    Parameters
    ----------
    l : int, optional
        l value for the first solid harmonic
    m : int, optional
        m value for the first solid harmonic
    lp : int, optional
        l' value for the second solid harmonic
    mp : int, optional
        m' value for the second solid harmonic
    t : int, optional
        t index for a given Hermite polynomial
    u : int, optional
        u index for a given Hermite polynomial
    v : int, optional
        v index for a given Hermite polynomial
    p : sympy symbol, optional
        symbolic representation of p, defined as alpha+beta where
        alpha is the exponent asigned to the first gtf and beta to
        the second gtf.
    Dx : sympy symbol, optional
        symbolic representation of Dx, defined as the x cartesian
        component for the difference between the coordinates for
        center A of the first gtf and the center P which results
        from the gaussian product rule.
    Dy : sympy symbol, optional
        symbolic representation of Dy, defined as the y cartesian
        component for the difference between the coordinates for
        center A of the first gtf and the center P which results
        from the gaussian product rule.
    Dz : sympy symbol, optional
        symbolic representation of Dz, defined as the z cartesian
        component for the difference between the coordinates for
        center A of the first gtf and the center P which results
        from the gaussian product rule.
    Dr2 : sympy symbol, optional
        symbolic representation of Dr^2, defined as the distance
        squared between the center A of the first gtf and the center
        P which results from the gaussian product rule.
    Dxp : sympy symbol, optional
        symbolic representation of Dx', defined as the x cartesian
        component for the difference between the coordinates for
        center B of the second gtf and the center P which results
        from the gaussian product rule.
    Dyp : sympy symbol, optional
        symbolic representation of Dy', defined as the y cartesian
        component for the difference between the coordinates for
        center B of the second gtf and the center P which results
        from the gaussian product rule.
    Dzp : sympy symbol, optional
        symbolic representation of Dz', defined as the z cartesian
        component for the difference between the coordinates for
        center B of the second gtf and the center P which results
        from the gaussian product rule.
    Drp2 : sympy symbol, optional
        symbolic representation of Dr'^2, defined as the distance
        squared between the center B of the second gtf and the center
        P which results from the gaussian product rule.

    """
    if(t < 0 or u < 0 or v < 0 or t + u + v > l + lp
       or l < abs(m) or lp < abs(mp)):
       #Invalid combinations of indices return zero
        return 0
    elif(t == 0 and u == 0 and v == 0 and l == 0 and lp == 0):
        #Starting coefficient for s-s product
        return 1
    elif(l >= lp):
        #When l is bigger than l', use recurrence relationships that
        #affects the value of l.
        if(l == abs(m)):
            if(l>1):
                if(m > 0):
                    return ((2*(l-1)+1)
                        * ((ECoeff(l-1,l-1,lp,mp,t-1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        - (ECoeff(l-1,-l+1,lp,mp,t,u-1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        + Dx * ECoeff(l-1,l-1,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2) - Dy
                        * ECoeff(l-1,-l+1,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        + (t+1)*ECoeff(l-1,l-1,lp,mp,t+1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        - (u+1)*ECoeff(l-1,-l+1,lp,mp,t,u+1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
                else:
                    return ((2*(l-1)+1)
                        * ((ECoeff(l-1,-l+1,lp,mp,t-1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        + (ECoeff(l-1,l-1,lp,mp,t,u-1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        + Dx * ECoeff(l-1,-l+1,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2) + Dy
                        * ECoeff(l-1,l-1,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        + (t+1)*ECoeff(l-1,-l+1,lp,mp,t+1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        + (u+1)*ECoeff(l-1,l-1,lp,mp,t,u+1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
            else:
               if(m > 0):
                    return ((2*(l-1)+1)
                       * ((ECoeff(l-1,l-1,lp,mp,t-1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                       + Dx * ECoeff(l-1,l-1,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                       + (t+1)*ECoeff(l-1,l-1,lp,mp,t+1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
               else:
                   return ((2*(l-1)+1)
                       * ((ECoeff(l-1,l-1,lp,mp,t,u-1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                       + Dy * ECoeff(l-1,l-1,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                       + (u+1)*ECoeff(l-1,l-1,lp,mp,t,u+1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
        else:
            return ((1/((l - 1) - abs(m) + 1)) * ((2*(l - 1) + 1)
                * ((ECoeff(l - 1, m, lp, mp, t, u, v - 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p)) + (Dz
                * ECoeff(l - 1, m, lp, mp, t, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) + ((v + 1)
                * ECoeff(l - 1, m, lp, mp, t, u, v + 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))) - ((l - 1
                + abs(m))*((ECoeff(l - 2, m, lp, mp, t - 2, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + ECoeff(l - 2, m, lp, mp, t, u - 2, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + ECoeff(l - 2, m, lp, mp, t, u, v - 2,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))/((2*p)**2)
                + ((Dx*ECoeff(l - 2, m, lp, mp, t - 1, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + Dy*ECoeff(l - 2, m, lp, mp, t, u - 1, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + Dz*ECoeff(l - 2, m, lp, mp, t, u, v - 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))/p)
                + ((Dr2 + (2*t + 2*u + 2*v + 3)/(2*p))
                * ECoeff(l - 2, m, lp, mp, t, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) + (2*Dx*(t + 1)
                * ECoeff(l - 2, m, lp, mp, t + 1, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) + (2*Dy*(u + 1)
                * ECoeff(l - 2, m, lp, mp, t, u + 1, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) + (2*Dz*(v + 1)
                * ECoeff(l - 2, m, lp, mp, t, u, v + 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) +((t + 2)*(t + 1)
                * ECoeff(l - 2, m, lp, mp, t + 2, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) + ((u + 2)*(u + 1)
                * ECoeff(l - 2, m, lp, mp, t, u + 2, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)) + ((v + 2)*(v + 1)
                * ECoeff(l - 2, m, lp, mp, t, u, v + 2,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))))))
    else:
        if(lp == abs(mp)):
            if(lp>1):
                if(mp>0):
                    return ((2*(lp-1)+1)
                        * ((ECoeff(l,m,lp-1,lp-1,t-1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        -(ECoeff(l,m,lp-1,-lp+1,t,u-1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        + Dxp * ECoeff(l,m,lp-1,lp-1,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2) - Dyp
                        * ECoeff(l,m,lp-1,-lp+1,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        + (t+1)*ECoeff(l,m,lp-1,lp-1,t+1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        - (u+1)*ECoeff(l,m,lp-1,-lp+1,t,u+1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
                else:
                    return ((2*(lp-1)+1)
                        * ((ECoeff(l,m,lp-1,-lp+1,t-1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        +(ECoeff(l,m,lp-1,lp-1,t,u-1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                        + Dxp * ECoeff(l,m,lp-1,-lp+1,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2) + Dyp
                        * ECoeff(l,m,lp-1,lp-1,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        + (t+1)*ECoeff(l,m,lp-1,-lp+1,t+1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                        + (u+1)*ECoeff(l,m,lp-1,lp-1,t,u+1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
            else:
                if(mp > 0):
                    return ((2*(lp-1)+1)
                       * ((ECoeff(l,m,lp-1,lp-1,t-1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                       + Dxp * ECoeff(l,m,lp-1,lp-1,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                       + (t+1)*ECoeff(l,m,lp-1,lp-1,t+1,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
                else:
                    return ((2*(lp-1)+1)
                       * ((ECoeff(l,m,lp-1,lp-1,t,u-1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                       + Dyp * ECoeff(l,m,lp-1,lp-1,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                       + (u+1)*ECoeff(l,m,lp-1,lp-1,t,u+1,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))

        else:
            return ((1/((lp - 1) - abs(mp) + 1)) * ((2*(lp - 1) + 1)
                *((ECoeff(l,m,lp - 1, mp, t, u, v - 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)/(2*p))
                + (Dzp * ECoeff(l,m,lp - 1, mp, t, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + ((v + 1) * ECoeff(l,m,lp - 1, mp, t, u, v + 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
                - ((lp - 1 + abs(mp))
                * ((ECoeff(l,m,lp - 2, mp, t - 2, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + ECoeff(l,m,lp - 2, mp, t, u - 2, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + ECoeff(l,m,lp - 2, mp, t, u, v - 2,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))/((2*p)**2)
                + ((Dxp*ECoeff(l,m,lp - 2, mp, t - 1, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + Dyp*ECoeff(l,m,lp - 2, mp, t, u - 1, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)
                + Dzp*ECoeff(l,m,lp - 2, mp, t, u, v - 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))/p)
                + ((Drp2 + (2*t + 2*u + 2*v + 3)/(2*p))
                * ECoeff(l,m,lp - 2, mp, t, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + (2*Dxp*(t + 1) * ECoeff(l,m,lp - 2, mp, t + 1, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + (2*Dyp*(u + 1) * ECoeff(l,m,lp - 2, mp, t, u + 1, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + (2*Dzp*(v + 1) * ECoeff(l,m,lp - 2, mp, t, u, v + 1,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + ((t + 2)*(t + 1) * ECoeff(l,m,lp - 2, mp, t + 2, u, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + ((u + 2)*(u + 1) * ECoeff(l,m,lp - 2, mp, t, u + 2, v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))
                + ((v + 2)*(v + 1) * ECoeff(l,m,lp - 2, mp, t, u, v + 2,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2))))))
