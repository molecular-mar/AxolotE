'''
A routine for generating the C++ code for evaluating the E coefficients
is here presented.
'''
from ERecurrence import *
from sympy import simplify
from sympy.printing import ccode
import shutil

# Create new file from header starter file
shutil.copyfile('starters/ECoeff_starter.cpp', 'AxolotE_Py.cpp')

# Variables names
p = symbols("expGamma")
Dx = symbols("xPtoA")
Dy = symbols("yPtoA")
Dz = symbols("zPtoA")
Dr2 = symbols("d2PtoA")
Dxp = symbols("xPtoB")
Dyp = symbols("yPtoB")
Dzp = symbols("zPtoB")
Drp2 = symbols("d2PtoB")

# Set maximum l values to generate
l_max = 2
lp_max = 2

pragma_line = "#pragma acc routine seq"
base_name = "ECoeff"

function_names = []
with open("ECoeffPy_Generated.cpp", 'a') as genCode:
    for l in range(l_max+1):
        for lp in range(lp_max+1):
            nCoeff = 0
            function_name = ("\nstd::vector<double> " + base_name + str(l) + "And" + str(lp)
                + "(double expAlpha, double expBeta,\n\tstd::array<double,3> \
muCoords, std::array<double,3> nuCoords)")
            function_names.append(function_name)
            genCode.write("\n" + pragma_line
                + "\nstd::vector<double> " + base_name + str(l) + "And" + str(lp)
                + "(double expAlpha, double expBeta,\n\tstd::array<double,3> \
muCoords, std::array<double,3> nuCoords)")
            genCode.write("\n{\
\n\tdouble expGamma{expAlpha + expBeta};\
\n\tdouble betaOverGamma{expBeta / expGamma};\
\n\tdouble xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};\
\n\tdouble yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};\
\n\tdouble zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};\
\n\tdouble d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};\
\n\tdouble alphaOverGamma{-expAlpha / expGamma};\
\n\tdouble xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};\
\n\tdouble yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};\
\n\tdouble zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};\
\n\tdouble d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};\
\n\tint nCoeffs {" + str(int((2*l + 1)*(2*lp + 1)*(l + lp + 1)*(l + lp +
   2)*(l + lp + 3)/6)) + "};\
\n\tstd::vector<double> expCoeffs(nCoeffs,0);")
            for m in range(-l, l+1):
                for mp in range(-lp, lp+1):
                    for t in range(l+lp+1):
                        for u in range(l+lp+1):
                            for v in range(l+lp+1):
                                if (t + u + v) <= (l + lp):
                                    m_strings = str(m) + str(mp)
                                    tuv_strings = str(t) + str(u) + str(v)
                                    genCode.write("\n\t// m,mp=" + m_strings + " t,u,v=" + tuv_strings)
                                    genCode.write("\n\texpCoeffs[" + str(nCoeff) + "] = "
                                            + ccode(simplify(ECoeff(l,m,lp,mp,t,u,v,p,Dx,Dy,Dz,Dr2,Dxp,Dyp,Dzp,Drp2)))
                                            +';')
                                    nCoeff += 1
            genCode.write("\n\treturn expCoeffs;\n}\n")

# Header files generation
shutil.copyfile('starters/ECoeff_starter.h', 'AxolotE_Py.h')
with open("AxolotePy.h", 'a') as genHeader:
    for function_name in function_names:
        genHeader.write(function_name + ";")
    genHeader.write("\n#endif")

