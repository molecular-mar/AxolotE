'''
This file prints the algebraic expressions generated using
Mathematica format, which can be used for their further
processing inside Mathematica tools.
'''

from ERecurrence import *
from sympy import simplify, mathematica_code as mcode

# Generate the expressions from l,lp=00 (ss) up to a given
# l_max and lp_max value.
l_max = 4
lp_max = 4

for l in range(l_max+1):
    for lp in range(lp_max+1):
        for m in range(-l, l+1):
            for mp in range(-lp, lp+1):
                for v in range(l+lp+1):
                    for u in range(l+lp+1):
                        for t in range(l+lp+1):
                            if (t + u + v) <= (l + lp):
                                print(mcode(simplify(ECoeff(l,m,lp,mp,t,u,v))))
