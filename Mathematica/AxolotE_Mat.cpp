#include "AxolotE_Mat.h"
#include <iostream>
#include <vector>
#include <array>
#include <cmath>

/*Functions for the generation of E coefficients used in the
expansion of the product of two Real Spherical Harmonics GTOs
into a single Hermite GTO. The functions here stored were
generated using a Mathematica notebook.*/


#pragma acc routine seq
std::vector<double> ECoeff0And0(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {1};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=00 t,u,v=000
	expCoeffs[0] = 1;
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff0And1(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {12};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=0-1 t,u,v=000
	expCoeffs[0] = yPtoB;
	// m,mp=0-1 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=0-1 t,u,v=010
	expCoeffs[2] = 1/(2.*expGamma);
	// m,mp=0-1 t,u,v=001
	expCoeffs[3] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[4] = zPtoB;
	// m,mp=00 t,u,v=100
	expCoeffs[5] = 1/(2.*expGamma);
	// m,mp=00 t,u,v=010
	expCoeffs[6] = 0;
	// m,mp=00 t,u,v=001
	expCoeffs[7] = 0;
	// m,mp=01 t,u,v=000
	expCoeffs[8] = xPtoB;
	// m,mp=01 t,u,v=100
	expCoeffs[9] = 0;
	// m,mp=01 t,u,v=010
	expCoeffs[10] = 0;
	// m,mp=01 t,u,v=001
	expCoeffs[11] = 1/(2.*expGamma);
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff0And2(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {50};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=0-2 t,u,v=000
	expCoeffs[0] = 6*xPtoB*yPtoB;
	// m,mp=0-2 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=0-2 t,u,v=200
	expCoeffs[2] = 0;
	// m,mp=0-2 t,u,v=010
	expCoeffs[3] = (3*xPtoB)/expGamma;
	// m,mp=0-2 t,u,v=110
	expCoeffs[4] = 0;
	// m,mp=0-2 t,u,v=020
	expCoeffs[5] = 0;
	// m,mp=0-2 t,u,v=001
	expCoeffs[6] = (3*yPtoB)/expGamma;
	// m,mp=0-2 t,u,v=101
	expCoeffs[7] = 0;
	// m,mp=0-2 t,u,v=011
	expCoeffs[8] = 3/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=002
	expCoeffs[9] = 0;
	// m,mp=0-1 t,u,v=000
	expCoeffs[10] = 3*yPtoB*zPtoB;
	// m,mp=0-1 t,u,v=100
	expCoeffs[11] = (3*yPtoB)/(2.*expGamma);
	// m,mp=0-1 t,u,v=200
	expCoeffs[12] = 0;
	// m,mp=0-1 t,u,v=010
	expCoeffs[13] = (3*zPtoB)/(2.*expGamma);
	// m,mp=0-1 t,u,v=110
	expCoeffs[14] = 3/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=020
	expCoeffs[15] = 0;
	// m,mp=0-1 t,u,v=001
	expCoeffs[16] = 0;
	// m,mp=0-1 t,u,v=101
	expCoeffs[17] = 0;
	// m,mp=0-1 t,u,v=011
	expCoeffs[18] = 0;
	// m,mp=0-1 t,u,v=002
	expCoeffs[19] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[20] = (-d2PtoB + 3*std::pow(zPtoB,2))/2.;
	// m,mp=00 t,u,v=100
	expCoeffs[21] = zPtoB/expGamma;
	// m,mp=00 t,u,v=200
	expCoeffs[22] = 1/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=010
	expCoeffs[23] = -0.5*yPtoB/expGamma;
	// m,mp=00 t,u,v=110
	expCoeffs[24] = 0;
	// m,mp=00 t,u,v=020
	expCoeffs[25] = -0.125*1/std::pow(expGamma,2);
	// m,mp=00 t,u,v=001
	expCoeffs[26] = -0.5*xPtoB/expGamma;
	// m,mp=00 t,u,v=101
	expCoeffs[27] = 0;
	// m,mp=00 t,u,v=011
	expCoeffs[28] = 0;
	// m,mp=00 t,u,v=002
	expCoeffs[29] = -0.125*1/std::pow(expGamma,2);
	// m,mp=01 t,u,v=000
	expCoeffs[30] = 3*xPtoB*zPtoB;
	// m,mp=01 t,u,v=100
	expCoeffs[31] = (3*xPtoB)/(2.*expGamma);
	// m,mp=01 t,u,v=200
	expCoeffs[32] = 0;
	// m,mp=01 t,u,v=010
	expCoeffs[33] = 0;
	// m,mp=01 t,u,v=110
	expCoeffs[34] = 0;
	// m,mp=01 t,u,v=020
	expCoeffs[35] = 0;
	// m,mp=01 t,u,v=001
	expCoeffs[36] = (3*zPtoB)/(2.*expGamma);
	// m,mp=01 t,u,v=101
	expCoeffs[37] = 3/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=011
	expCoeffs[38] = 0;
	// m,mp=01 t,u,v=002
	expCoeffs[39] = 0;
	// m,mp=02 t,u,v=000
	expCoeffs[40] = 3*(xPtoB - yPtoB)*(xPtoB + yPtoB);
	// m,mp=02 t,u,v=100
	expCoeffs[41] = 0;
	// m,mp=02 t,u,v=200
	expCoeffs[42] = 0;
	// m,mp=02 t,u,v=010
	expCoeffs[43] = (-3*yPtoB)/expGamma;
	// m,mp=02 t,u,v=110
	expCoeffs[44] = 0;
	// m,mp=02 t,u,v=020
	expCoeffs[45] = -3/(4.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=001
	expCoeffs[46] = (3*xPtoB)/expGamma;
	// m,mp=02 t,u,v=101
	expCoeffs[47] = 0;
	// m,mp=02 t,u,v=011
	expCoeffs[48] = 0;
	// m,mp=02 t,u,v=002
	expCoeffs[49] = 3/(4.*std::pow(expGamma,2));
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff1And0(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {12};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=-10 t,u,v=000
	expCoeffs[0] = yPtoA;
	// m,mp=-10 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=-10 t,u,v=010
	expCoeffs[2] = 1/(2.*expGamma);
	// m,mp=-10 t,u,v=001
	expCoeffs[3] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[4] = zPtoA;
	// m,mp=00 t,u,v=100
	expCoeffs[5] = 1/(2.*expGamma);
	// m,mp=00 t,u,v=010
	expCoeffs[6] = 0;
	// m,mp=00 t,u,v=001
	expCoeffs[7] = 0;
	// m,mp=10 t,u,v=000
	expCoeffs[8] = xPtoA;
	// m,mp=10 t,u,v=100
	expCoeffs[9] = 0;
	// m,mp=10 t,u,v=010
	expCoeffs[10] = 0;
	// m,mp=10 t,u,v=001
	expCoeffs[11] = 1/(2.*expGamma);
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff1And1(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {90};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=-1-1 t,u,v=000
	expCoeffs[0] = 1/(2.*expGamma) + yPtoA*yPtoB;
	// m,mp=-1-1 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=-1-1 t,u,v=200
	expCoeffs[2] = 0;
	// m,mp=-1-1 t,u,v=010
	expCoeffs[3] = (yPtoA + yPtoB)/(2.*expGamma);
	// m,mp=-1-1 t,u,v=110
	expCoeffs[4] = 0;
	// m,mp=-1-1 t,u,v=020
	expCoeffs[5] = 1/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=001
	expCoeffs[6] = 0;
	// m,mp=-1-1 t,u,v=101
	expCoeffs[7] = 0;
	// m,mp=-1-1 t,u,v=011
	expCoeffs[8] = 0;
	// m,mp=-1-1 t,u,v=002
	expCoeffs[9] = 0;
	// m,mp=-10 t,u,v=000
	expCoeffs[10] = yPtoA*zPtoB;
	// m,mp=-10 t,u,v=100
	expCoeffs[11] = yPtoA/(2.*expGamma);
	// m,mp=-10 t,u,v=200
	expCoeffs[12] = 0;
	// m,mp=-10 t,u,v=010
	expCoeffs[13] = zPtoB/(2.*expGamma);
	// m,mp=-10 t,u,v=110
	expCoeffs[14] = 1/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=020
	expCoeffs[15] = 0;
	// m,mp=-10 t,u,v=001
	expCoeffs[16] = 0;
	// m,mp=-10 t,u,v=101
	expCoeffs[17] = 0;
	// m,mp=-10 t,u,v=011
	expCoeffs[18] = 0;
	// m,mp=-10 t,u,v=002
	expCoeffs[19] = 0;
	// m,mp=-11 t,u,v=000
	expCoeffs[20] = xPtoB*yPtoA;
	// m,mp=-11 t,u,v=100
	expCoeffs[21] = 0;
	// m,mp=-11 t,u,v=200
	expCoeffs[22] = 0;
	// m,mp=-11 t,u,v=010
	expCoeffs[23] = xPtoB/(2.*expGamma);
	// m,mp=-11 t,u,v=110
	expCoeffs[24] = 0;
	// m,mp=-11 t,u,v=020
	expCoeffs[25] = 0;
	// m,mp=-11 t,u,v=001
	expCoeffs[26] = yPtoA/(2.*expGamma);
	// m,mp=-11 t,u,v=101
	expCoeffs[27] = 0;
	// m,mp=-11 t,u,v=011
	expCoeffs[28] = 1/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=002
	expCoeffs[29] = 0;
	// m,mp=0-1 t,u,v=000
	expCoeffs[30] = yPtoB*zPtoA;
	// m,mp=0-1 t,u,v=100
	expCoeffs[31] = yPtoB/(2.*expGamma);
	// m,mp=0-1 t,u,v=200
	expCoeffs[32] = 0;
	// m,mp=0-1 t,u,v=010
	expCoeffs[33] = zPtoA/(2.*expGamma);
	// m,mp=0-1 t,u,v=110
	expCoeffs[34] = 1/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=020
	expCoeffs[35] = 0;
	// m,mp=0-1 t,u,v=001
	expCoeffs[36] = 0;
	// m,mp=0-1 t,u,v=101
	expCoeffs[37] = 0;
	// m,mp=0-1 t,u,v=011
	expCoeffs[38] = 0;
	// m,mp=0-1 t,u,v=002
	expCoeffs[39] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[40] = 1/(2.*expGamma) + zPtoA*zPtoB;
	// m,mp=00 t,u,v=100
	expCoeffs[41] = (zPtoA + zPtoB)/(2.*expGamma);
	// m,mp=00 t,u,v=200
	expCoeffs[42] = 1/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=010
	expCoeffs[43] = 0;
	// m,mp=00 t,u,v=110
	expCoeffs[44] = 0;
	// m,mp=00 t,u,v=020
	expCoeffs[45] = 0;
	// m,mp=00 t,u,v=001
	expCoeffs[46] = 0;
	// m,mp=00 t,u,v=101
	expCoeffs[47] = 0;
	// m,mp=00 t,u,v=011
	expCoeffs[48] = 0;
	// m,mp=00 t,u,v=002
	expCoeffs[49] = 0;
	// m,mp=01 t,u,v=000
	expCoeffs[50] = xPtoB*zPtoA;
	// m,mp=01 t,u,v=100
	expCoeffs[51] = xPtoB/(2.*expGamma);
	// m,mp=01 t,u,v=200
	expCoeffs[52] = 0;
	// m,mp=01 t,u,v=010
	expCoeffs[53] = 0;
	// m,mp=01 t,u,v=110
	expCoeffs[54] = 0;
	// m,mp=01 t,u,v=020
	expCoeffs[55] = 0;
	// m,mp=01 t,u,v=001
	expCoeffs[56] = zPtoA/(2.*expGamma);
	// m,mp=01 t,u,v=101
	expCoeffs[57] = 1/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=011
	expCoeffs[58] = 0;
	// m,mp=01 t,u,v=002
	expCoeffs[59] = 0;
	// m,mp=1-1 t,u,v=000
	expCoeffs[60] = xPtoA*yPtoB;
	// m,mp=1-1 t,u,v=100
	expCoeffs[61] = 0;
	// m,mp=1-1 t,u,v=200
	expCoeffs[62] = 0;
	// m,mp=1-1 t,u,v=010
	expCoeffs[63] = xPtoA/(2.*expGamma);
	// m,mp=1-1 t,u,v=110
	expCoeffs[64] = 0;
	// m,mp=1-1 t,u,v=020
	expCoeffs[65] = 0;
	// m,mp=1-1 t,u,v=001
	expCoeffs[66] = yPtoB/(2.*expGamma);
	// m,mp=1-1 t,u,v=101
	expCoeffs[67] = 0;
	// m,mp=1-1 t,u,v=011
	expCoeffs[68] = 1/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=002
	expCoeffs[69] = 0;
	// m,mp=10 t,u,v=000
	expCoeffs[70] = xPtoA*zPtoB;
	// m,mp=10 t,u,v=100
	expCoeffs[71] = xPtoA/(2.*expGamma);
	// m,mp=10 t,u,v=200
	expCoeffs[72] = 0;
	// m,mp=10 t,u,v=010
	expCoeffs[73] = 0;
	// m,mp=10 t,u,v=110
	expCoeffs[74] = 0;
	// m,mp=10 t,u,v=020
	expCoeffs[75] = 0;
	// m,mp=10 t,u,v=001
	expCoeffs[76] = zPtoB/(2.*expGamma);
	// m,mp=10 t,u,v=101
	expCoeffs[77] = 1/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=011
	expCoeffs[78] = 0;
	// m,mp=10 t,u,v=002
	expCoeffs[79] = 0;
	// m,mp=11 t,u,v=000
	expCoeffs[80] = 1/(2.*expGamma) + xPtoA*xPtoB;
	// m,mp=11 t,u,v=100
	expCoeffs[81] = 0;
	// m,mp=11 t,u,v=200
	expCoeffs[82] = 0;
	// m,mp=11 t,u,v=010
	expCoeffs[83] = 0;
	// m,mp=11 t,u,v=110
	expCoeffs[84] = 0;
	// m,mp=11 t,u,v=020
	expCoeffs[85] = 0;
	// m,mp=11 t,u,v=001
	expCoeffs[86] = (xPtoA + xPtoB)/(2.*expGamma);
	// m,mp=11 t,u,v=101
	expCoeffs[87] = 0;
	// m,mp=11 t,u,v=011
	expCoeffs[88] = 0;
	// m,mp=11 t,u,v=002
	expCoeffs[89] = 1/(4.*std::pow(expGamma,2));
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff1And2(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {300};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=-1-2 t,u,v=000
	expCoeffs[0] = 3*xPtoB*(1/expGamma + 2*yPtoA*yPtoB);
	// m,mp=-1-2 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=-1-2 t,u,v=200
	expCoeffs[2] = 0;
	// m,mp=-1-2 t,u,v=300
	expCoeffs[3] = 0;
	// m,mp=-1-2 t,u,v=010
	expCoeffs[4] = (3*xPtoB*(yPtoA + yPtoB))/expGamma;
	// m,mp=-1-2 t,u,v=110
	expCoeffs[5] = 0;
	// m,mp=-1-2 t,u,v=210
	expCoeffs[6] = 0;
	// m,mp=-1-2 t,u,v=020
	expCoeffs[7] = (3*xPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=120
	expCoeffs[8] = 0;
	// m,mp=-1-2 t,u,v=030
	expCoeffs[9] = 0;
	// m,mp=-1-2 t,u,v=001
	expCoeffs[10] = (3 + 6*expGamma*yPtoA*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=101
	expCoeffs[11] = 0;
	// m,mp=-1-2 t,u,v=201
	expCoeffs[12] = 0;
	// m,mp=-1-2 t,u,v=011
	expCoeffs[13] = (3*(yPtoA + yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=111
	expCoeffs[14] = 0;
	// m,mp=-1-2 t,u,v=021
	expCoeffs[15] = 3/(4.*std::pow(expGamma,3));
	// m,mp=-1-2 t,u,v=002
	expCoeffs[16] = 0;
	// m,mp=-1-2 t,u,v=102
	expCoeffs[17] = 0;
	// m,mp=-1-2 t,u,v=012
	expCoeffs[18] = 0;
	// m,mp=-1-2 t,u,v=003
	expCoeffs[19] = 0;
	// m,mp=-1-1 t,u,v=000
	expCoeffs[20] = (3*(1/expGamma + 2*yPtoA*yPtoB)*zPtoB)/2.;
	// m,mp=-1-1 t,u,v=100
	expCoeffs[21] = (3 + 6*expGamma*yPtoA*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=200
	expCoeffs[22] = 0;
	// m,mp=-1-1 t,u,v=300
	expCoeffs[23] = 0;
	// m,mp=-1-1 t,u,v=010
	expCoeffs[24] = (3*(yPtoA + yPtoB)*zPtoB)/(2.*expGamma);
	// m,mp=-1-1 t,u,v=110
	expCoeffs[25] = (3*(yPtoA + yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=210
	expCoeffs[26] = 0;
	// m,mp=-1-1 t,u,v=020
	expCoeffs[27] = (3*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=120
	expCoeffs[28] = 3/(8.*std::pow(expGamma,3));
	// m,mp=-1-1 t,u,v=030
	expCoeffs[29] = 0;
	// m,mp=-1-1 t,u,v=001
	expCoeffs[30] = 0;
	// m,mp=-1-1 t,u,v=101
	expCoeffs[31] = 0;
	// m,mp=-1-1 t,u,v=201
	expCoeffs[32] = 0;
	// m,mp=-1-1 t,u,v=011
	expCoeffs[33] = 0;
	// m,mp=-1-1 t,u,v=111
	expCoeffs[34] = 0;
	// m,mp=-1-1 t,u,v=021
	expCoeffs[35] = 0;
	// m,mp=-1-1 t,u,v=002
	expCoeffs[36] = 0;
	// m,mp=-1-1 t,u,v=102
	expCoeffs[37] = 0;
	// m,mp=-1-1 t,u,v=012
	expCoeffs[38] = 0;
	// m,mp=-1-1 t,u,v=003
	expCoeffs[39] = 0;
	// m,mp=-10 t,u,v=000
	expCoeffs[40] = -0.5*(yPtoB + expGamma*yPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))/expGamma;
	// m,mp=-10 t,u,v=100
	expCoeffs[41] = (yPtoA*zPtoB)/expGamma;
	// m,mp=-10 t,u,v=200
	expCoeffs[42] = yPtoA/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=300
	expCoeffs[43] = 0;
	// m,mp=-10 t,u,v=010
	expCoeffs[44] = -0.25*(1 + expGamma*(d2PtoB + 2*yPtoA*yPtoB - 3*std::pow(zPtoB,2)))/std::pow(expGamma,2);
	// m,mp=-10 t,u,v=110
	expCoeffs[45] = zPtoB/(2.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=210
	expCoeffs[46] = 1/(8.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=020
	expCoeffs[47] = -0.125*(yPtoA + 2*yPtoB)/std::pow(expGamma,2);
	// m,mp=-10 t,u,v=120
	expCoeffs[48] = 0;
	// m,mp=-10 t,u,v=030
	expCoeffs[49] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=-10 t,u,v=001
	expCoeffs[50] = -0.5*(xPtoB*yPtoA)/expGamma;
	// m,mp=-10 t,u,v=101
	expCoeffs[51] = 0;
	// m,mp=-10 t,u,v=201
	expCoeffs[52] = 0;
	// m,mp=-10 t,u,v=011
	expCoeffs[53] = -0.25*xPtoB/std::pow(expGamma,2);
	// m,mp=-10 t,u,v=111
	expCoeffs[54] = 0;
	// m,mp=-10 t,u,v=021
	expCoeffs[55] = 0;
	// m,mp=-10 t,u,v=002
	expCoeffs[56] = -0.125*yPtoA/std::pow(expGamma,2);
	// m,mp=-10 t,u,v=102
	expCoeffs[57] = 0;
	// m,mp=-10 t,u,v=012
	expCoeffs[58] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=-10 t,u,v=003
	expCoeffs[59] = 0;
	// m,mp=-11 t,u,v=000
	expCoeffs[60] = 3*xPtoB*yPtoA*zPtoB;
	// m,mp=-11 t,u,v=100
	expCoeffs[61] = (3*xPtoB*yPtoA)/(2.*expGamma);
	// m,mp=-11 t,u,v=200
	expCoeffs[62] = 0;
	// m,mp=-11 t,u,v=300
	expCoeffs[63] = 0;
	// m,mp=-11 t,u,v=010
	expCoeffs[64] = (3*xPtoB*zPtoB)/(2.*expGamma);
	// m,mp=-11 t,u,v=110
	expCoeffs[65] = (3*xPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=210
	expCoeffs[66] = 0;
	// m,mp=-11 t,u,v=020
	expCoeffs[67] = 0;
	// m,mp=-11 t,u,v=120
	expCoeffs[68] = 0;
	// m,mp=-11 t,u,v=030
	expCoeffs[69] = 0;
	// m,mp=-11 t,u,v=001
	expCoeffs[70] = (3*yPtoA*zPtoB)/(2.*expGamma);
	// m,mp=-11 t,u,v=101
	expCoeffs[71] = (3*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=201
	expCoeffs[72] = 0;
	// m,mp=-11 t,u,v=011
	expCoeffs[73] = (3*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=111
	expCoeffs[74] = 3/(8.*std::pow(expGamma,3));
	// m,mp=-11 t,u,v=021
	expCoeffs[75] = 0;
	// m,mp=-11 t,u,v=002
	expCoeffs[76] = 0;
	// m,mp=-11 t,u,v=102
	expCoeffs[77] = 0;
	// m,mp=-11 t,u,v=012
	expCoeffs[78] = 0;
	// m,mp=-11 t,u,v=003
	expCoeffs[79] = 0;
	// m,mp=-12 t,u,v=000
	expCoeffs[80] = 3*std::pow(xPtoB,2)*yPtoA - (3*yPtoB*(1 + expGamma*yPtoA*yPtoB))/expGamma;
	// m,mp=-12 t,u,v=100
	expCoeffs[81] = 0;
	// m,mp=-12 t,u,v=200
	expCoeffs[82] = 0;
	// m,mp=-12 t,u,v=300
	expCoeffs[83] = 0;
	// m,mp=-12 t,u,v=010
	expCoeffs[84] = (-3 + 3*expGamma*std::pow(xPtoB,2) - 3*expGamma*yPtoB*(2*yPtoA + yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=110
	expCoeffs[85] = 0;
	// m,mp=-12 t,u,v=210
	expCoeffs[86] = 0;
	// m,mp=-12 t,u,v=020
	expCoeffs[87] = (-3*(yPtoA + 2*yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=120
	expCoeffs[88] = 0;
	// m,mp=-12 t,u,v=030
	expCoeffs[89] = -3/(8.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=001
	expCoeffs[90] = (3*xPtoB*yPtoA)/expGamma;
	// m,mp=-12 t,u,v=101
	expCoeffs[91] = 0;
	// m,mp=-12 t,u,v=201
	expCoeffs[92] = 0;
	// m,mp=-12 t,u,v=011
	expCoeffs[93] = (3*xPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=111
	expCoeffs[94] = 0;
	// m,mp=-12 t,u,v=021
	expCoeffs[95] = 0;
	// m,mp=-12 t,u,v=002
	expCoeffs[96] = (3*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=102
	expCoeffs[97] = 0;
	// m,mp=-12 t,u,v=012
	expCoeffs[98] = 3/(8.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=003
	expCoeffs[99] = 0;
	// m,mp=0-2 t,u,v=000
	expCoeffs[100] = 6*xPtoB*yPtoB*zPtoA;
	// m,mp=0-2 t,u,v=100
	expCoeffs[101] = (3*xPtoB*yPtoB)/expGamma;
	// m,mp=0-2 t,u,v=200
	expCoeffs[102] = 0;
	// m,mp=0-2 t,u,v=300
	expCoeffs[103] = 0;
	// m,mp=0-2 t,u,v=010
	expCoeffs[104] = (3*xPtoB*zPtoA)/expGamma;
	// m,mp=0-2 t,u,v=110
	expCoeffs[105] = (3*xPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=210
	expCoeffs[106] = 0;
	// m,mp=0-2 t,u,v=020
	expCoeffs[107] = 0;
	// m,mp=0-2 t,u,v=120
	expCoeffs[108] = 0;
	// m,mp=0-2 t,u,v=030
	expCoeffs[109] = 0;
	// m,mp=0-2 t,u,v=001
	expCoeffs[110] = (3*yPtoB*zPtoA)/expGamma;
	// m,mp=0-2 t,u,v=101
	expCoeffs[111] = (3*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=201
	expCoeffs[112] = 0;
	// m,mp=0-2 t,u,v=011
	expCoeffs[113] = (3*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=111
	expCoeffs[114] = 3/(4.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=021
	expCoeffs[115] = 0;
	// m,mp=0-2 t,u,v=002
	expCoeffs[116] = 0;
	// m,mp=0-2 t,u,v=102
	expCoeffs[117] = 0;
	// m,mp=0-2 t,u,v=012
	expCoeffs[118] = 0;
	// m,mp=0-2 t,u,v=003
	expCoeffs[119] = 0;
	// m,mp=0-1 t,u,v=000
	expCoeffs[120] = (3*yPtoB*(1/expGamma + 2*zPtoA*zPtoB))/2.;
	// m,mp=0-1 t,u,v=100
	expCoeffs[121] = (3*yPtoB*(zPtoA + zPtoB))/(2.*expGamma);
	// m,mp=0-1 t,u,v=200
	expCoeffs[122] = (3*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=300
	expCoeffs[123] = 0;
	// m,mp=0-1 t,u,v=010
	expCoeffs[124] = (3 + 6*expGamma*zPtoA*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=110
	expCoeffs[125] = (3*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=210
	expCoeffs[126] = 3/(8.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=020
	expCoeffs[127] = 0;
	// m,mp=0-1 t,u,v=120
	expCoeffs[128] = 0;
	// m,mp=0-1 t,u,v=030
	expCoeffs[129] = 0;
	// m,mp=0-1 t,u,v=001
	expCoeffs[130] = 0;
	// m,mp=0-1 t,u,v=101
	expCoeffs[131] = 0;
	// m,mp=0-1 t,u,v=201
	expCoeffs[132] = 0;
	// m,mp=0-1 t,u,v=011
	expCoeffs[133] = 0;
	// m,mp=0-1 t,u,v=111
	expCoeffs[134] = 0;
	// m,mp=0-1 t,u,v=021
	expCoeffs[135] = 0;
	// m,mp=0-1 t,u,v=002
	expCoeffs[136] = 0;
	// m,mp=0-1 t,u,v=102
	expCoeffs[137] = 0;
	// m,mp=0-1 t,u,v=012
	expCoeffs[138] = 0;
	// m,mp=0-1 t,u,v=003
	expCoeffs[139] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[140] = zPtoB/expGamma - (zPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))/2.;
	// m,mp=00 t,u,v=100
	expCoeffs[141] = (2 - d2PtoB*expGamma + expGamma*zPtoB*(4*zPtoA + 3*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=200
	expCoeffs[142] = (zPtoA + 2*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=300
	expCoeffs[143] = 1/(8.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=010
	expCoeffs[144] = -0.5*(yPtoB*zPtoA)/expGamma;
	// m,mp=00 t,u,v=110
	expCoeffs[145] = -0.25*yPtoB/std::pow(expGamma,2);
	// m,mp=00 t,u,v=210
	expCoeffs[146] = 0;
	// m,mp=00 t,u,v=020
	expCoeffs[147] = -0.125*zPtoA/std::pow(expGamma,2);
	// m,mp=00 t,u,v=120
	expCoeffs[148] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=00 t,u,v=030
	expCoeffs[149] = 0;
	// m,mp=00 t,u,v=001
	expCoeffs[150] = -0.5*(xPtoB*zPtoA)/expGamma;
	// m,mp=00 t,u,v=101
	expCoeffs[151] = -0.25*xPtoB/std::pow(expGamma,2);
	// m,mp=00 t,u,v=201
	expCoeffs[152] = 0;
	// m,mp=00 t,u,v=011
	expCoeffs[153] = 0;
	// m,mp=00 t,u,v=111
	expCoeffs[154] = 0;
	// m,mp=00 t,u,v=021
	expCoeffs[155] = 0;
	// m,mp=00 t,u,v=002
	expCoeffs[156] = -0.125*zPtoA/std::pow(expGamma,2);
	// m,mp=00 t,u,v=102
	expCoeffs[157] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=00 t,u,v=012
	expCoeffs[158] = 0;
	// m,mp=00 t,u,v=003
	expCoeffs[159] = 0;
	// m,mp=01 t,u,v=000
	expCoeffs[160] = (3*xPtoB*(1/expGamma + 2*zPtoA*zPtoB))/2.;
	// m,mp=01 t,u,v=100
	expCoeffs[161] = (3*xPtoB*(zPtoA + zPtoB))/(2.*expGamma);
	// m,mp=01 t,u,v=200
	expCoeffs[162] = (3*xPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=300
	expCoeffs[163] = 0;
	// m,mp=01 t,u,v=010
	expCoeffs[164] = 0;
	// m,mp=01 t,u,v=110
	expCoeffs[165] = 0;
	// m,mp=01 t,u,v=210
	expCoeffs[166] = 0;
	// m,mp=01 t,u,v=020
	expCoeffs[167] = 0;
	// m,mp=01 t,u,v=120
	expCoeffs[168] = 0;
	// m,mp=01 t,u,v=030
	expCoeffs[169] = 0;
	// m,mp=01 t,u,v=001
	expCoeffs[170] = (3 + 6*expGamma*zPtoA*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=101
	expCoeffs[171] = (3*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=201
	expCoeffs[172] = 3/(8.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=011
	expCoeffs[173] = 0;
	// m,mp=01 t,u,v=111
	expCoeffs[174] = 0;
	// m,mp=01 t,u,v=021
	expCoeffs[175] = 0;
	// m,mp=01 t,u,v=002
	expCoeffs[176] = 0;
	// m,mp=01 t,u,v=102
	expCoeffs[177] = 0;
	// m,mp=01 t,u,v=012
	expCoeffs[178] = 0;
	// m,mp=01 t,u,v=003
	expCoeffs[179] = 0;
	// m,mp=02 t,u,v=000
	expCoeffs[180] = 3*(xPtoB - yPtoB)*(xPtoB + yPtoB)*zPtoA;
	// m,mp=02 t,u,v=100
	expCoeffs[181] = (3*(xPtoB - yPtoB)*(xPtoB + yPtoB))/(2.*expGamma);
	// m,mp=02 t,u,v=200
	expCoeffs[182] = 0;
	// m,mp=02 t,u,v=300
	expCoeffs[183] = 0;
	// m,mp=02 t,u,v=010
	expCoeffs[184] = (-3*yPtoB*zPtoA)/expGamma;
	// m,mp=02 t,u,v=110
	expCoeffs[185] = (-3*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=210
	expCoeffs[186] = 0;
	// m,mp=02 t,u,v=020
	expCoeffs[187] = (-3*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=120
	expCoeffs[188] = -3/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=030
	expCoeffs[189] = 0;
	// m,mp=02 t,u,v=001
	expCoeffs[190] = (3*xPtoB*zPtoA)/expGamma;
	// m,mp=02 t,u,v=101
	expCoeffs[191] = (3*xPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=201
	expCoeffs[192] = 0;
	// m,mp=02 t,u,v=011
	expCoeffs[193] = 0;
	// m,mp=02 t,u,v=111
	expCoeffs[194] = 0;
	// m,mp=02 t,u,v=021
	expCoeffs[195] = 0;
	// m,mp=02 t,u,v=002
	expCoeffs[196] = (3*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=102
	expCoeffs[197] = 3/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=012
	expCoeffs[198] = 0;
	// m,mp=02 t,u,v=003
	expCoeffs[199] = 0;
	// m,mp=1-2 t,u,v=000
	expCoeffs[200] = 3*(1/expGamma + 2*xPtoA*xPtoB)*yPtoB;
	// m,mp=1-2 t,u,v=100
	expCoeffs[201] = 0;
	// m,mp=1-2 t,u,v=200
	expCoeffs[202] = 0;
	// m,mp=1-2 t,u,v=300
	expCoeffs[203] = 0;
	// m,mp=1-2 t,u,v=010
	expCoeffs[204] = (3 + 6*expGamma*xPtoA*xPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=110
	expCoeffs[205] = 0;
	// m,mp=1-2 t,u,v=210
	expCoeffs[206] = 0;
	// m,mp=1-2 t,u,v=020
	expCoeffs[207] = 0;
	// m,mp=1-2 t,u,v=120
	expCoeffs[208] = 0;
	// m,mp=1-2 t,u,v=030
	expCoeffs[209] = 0;
	// m,mp=1-2 t,u,v=001
	expCoeffs[210] = (3*(xPtoA + xPtoB)*yPtoB)/expGamma;
	// m,mp=1-2 t,u,v=101
	expCoeffs[211] = 0;
	// m,mp=1-2 t,u,v=201
	expCoeffs[212] = 0;
	// m,mp=1-2 t,u,v=011
	expCoeffs[213] = (3*(xPtoA + xPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=111
	expCoeffs[214] = 0;
	// m,mp=1-2 t,u,v=021
	expCoeffs[215] = 0;
	// m,mp=1-2 t,u,v=002
	expCoeffs[216] = (3*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=102
	expCoeffs[217] = 0;
	// m,mp=1-2 t,u,v=012
	expCoeffs[218] = 3/(4.*std::pow(expGamma,3));
	// m,mp=1-2 t,u,v=003
	expCoeffs[219] = 0;
	// m,mp=1-1 t,u,v=000
	expCoeffs[220] = 3*xPtoA*yPtoB*zPtoB;
	// m,mp=1-1 t,u,v=100
	expCoeffs[221] = (3*xPtoA*yPtoB)/(2.*expGamma);
	// m,mp=1-1 t,u,v=200
	expCoeffs[222] = 0;
	// m,mp=1-1 t,u,v=300
	expCoeffs[223] = 0;
	// m,mp=1-1 t,u,v=010
	expCoeffs[224] = (3*xPtoA*zPtoB)/(2.*expGamma);
	// m,mp=1-1 t,u,v=110
	expCoeffs[225] = (3*xPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=210
	expCoeffs[226] = 0;
	// m,mp=1-1 t,u,v=020
	expCoeffs[227] = 0;
	// m,mp=1-1 t,u,v=120
	expCoeffs[228] = 0;
	// m,mp=1-1 t,u,v=030
	expCoeffs[229] = 0;
	// m,mp=1-1 t,u,v=001
	expCoeffs[230] = (3*yPtoB*zPtoB)/(2.*expGamma);
	// m,mp=1-1 t,u,v=101
	expCoeffs[231] = (3*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=201
	expCoeffs[232] = 0;
	// m,mp=1-1 t,u,v=011
	expCoeffs[233] = (3*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=111
	expCoeffs[234] = 3/(8.*std::pow(expGamma,3));
	// m,mp=1-1 t,u,v=021
	expCoeffs[235] = 0;
	// m,mp=1-1 t,u,v=002
	expCoeffs[236] = 0;
	// m,mp=1-1 t,u,v=102
	expCoeffs[237] = 0;
	// m,mp=1-1 t,u,v=012
	expCoeffs[238] = 0;
	// m,mp=1-1 t,u,v=003
	expCoeffs[239] = 0;
	// m,mp=10 t,u,v=000
	expCoeffs[240] = -0.5*(xPtoB + expGamma*xPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))/expGamma;
	// m,mp=10 t,u,v=100
	expCoeffs[241] = (xPtoA*zPtoB)/expGamma;
	// m,mp=10 t,u,v=200
	expCoeffs[242] = xPtoA/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=300
	expCoeffs[243] = 0;
	// m,mp=10 t,u,v=010
	expCoeffs[244] = -0.5*(xPtoA*yPtoB)/expGamma;
	// m,mp=10 t,u,v=110
	expCoeffs[245] = 0;
	// m,mp=10 t,u,v=210
	expCoeffs[246] = 0;
	// m,mp=10 t,u,v=020
	expCoeffs[247] = -0.125*xPtoA/std::pow(expGamma,2);
	// m,mp=10 t,u,v=120
	expCoeffs[248] = 0;
	// m,mp=10 t,u,v=030
	expCoeffs[249] = 0;
	// m,mp=10 t,u,v=001
	expCoeffs[250] = -0.25*(1 + expGamma*(d2PtoB + 2*xPtoA*xPtoB - 3*std::pow(zPtoB,2)))/std::pow(expGamma,2);
	// m,mp=10 t,u,v=101
	expCoeffs[251] = zPtoB/(2.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=201
	expCoeffs[252] = 1/(8.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=011
	expCoeffs[253] = -0.25*yPtoB/std::pow(expGamma,2);
	// m,mp=10 t,u,v=111
	expCoeffs[254] = 0;
	// m,mp=10 t,u,v=021
	expCoeffs[255] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=10 t,u,v=002
	expCoeffs[256] = -0.125*(xPtoA + 2*xPtoB)/std::pow(expGamma,2);
	// m,mp=10 t,u,v=102
	expCoeffs[257] = 0;
	// m,mp=10 t,u,v=012
	expCoeffs[258] = 0;
	// m,mp=10 t,u,v=003
	expCoeffs[259] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=11 t,u,v=000
	expCoeffs[260] = (3*(1/expGamma + 2*xPtoA*xPtoB)*zPtoB)/2.;
	// m,mp=11 t,u,v=100
	expCoeffs[261] = (3 + 6*expGamma*xPtoA*xPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=200
	expCoeffs[262] = 0;
	// m,mp=11 t,u,v=300
	expCoeffs[263] = 0;
	// m,mp=11 t,u,v=010
	expCoeffs[264] = 0;
	// m,mp=11 t,u,v=110
	expCoeffs[265] = 0;
	// m,mp=11 t,u,v=210
	expCoeffs[266] = 0;
	// m,mp=11 t,u,v=020
	expCoeffs[267] = 0;
	// m,mp=11 t,u,v=120
	expCoeffs[268] = 0;
	// m,mp=11 t,u,v=030
	expCoeffs[269] = 0;
	// m,mp=11 t,u,v=001
	expCoeffs[270] = (3*(xPtoA + xPtoB)*zPtoB)/(2.*expGamma);
	// m,mp=11 t,u,v=101
	expCoeffs[271] = (3*(xPtoA + xPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=201
	expCoeffs[272] = 0;
	// m,mp=11 t,u,v=011
	expCoeffs[273] = 0;
	// m,mp=11 t,u,v=111
	expCoeffs[274] = 0;
	// m,mp=11 t,u,v=021
	expCoeffs[275] = 0;
	// m,mp=11 t,u,v=002
	expCoeffs[276] = (3*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=102
	expCoeffs[277] = 3/(8.*std::pow(expGamma,3));
	// m,mp=11 t,u,v=012
	expCoeffs[278] = 0;
	// m,mp=11 t,u,v=003
	expCoeffs[279] = 0;
	// m,mp=12 t,u,v=000
	expCoeffs[280] = (3*(xPtoB + expGamma*xPtoA*(xPtoB - yPtoB)*(xPtoB + yPtoB)))/expGamma;
	// m,mp=12 t,u,v=100
	expCoeffs[281] = 0;
	// m,mp=12 t,u,v=200
	expCoeffs[282] = 0;
	// m,mp=12 t,u,v=300
	expCoeffs[283] = 0;
	// m,mp=12 t,u,v=010
	expCoeffs[284] = (-3*xPtoA*yPtoB)/expGamma;
	// m,mp=12 t,u,v=110
	expCoeffs[285] = 0;
	// m,mp=12 t,u,v=210
	expCoeffs[286] = 0;
	// m,mp=12 t,u,v=020
	expCoeffs[287] = (-3*xPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=120
	expCoeffs[288] = 0;
	// m,mp=12 t,u,v=030
	expCoeffs[289] = 0;
	// m,mp=12 t,u,v=001
	expCoeffs[290] = (3*(1 + 2*expGamma*xPtoA*xPtoB + expGamma*(xPtoB - yPtoB)*(xPtoB + yPtoB)))/(2.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=101
	expCoeffs[291] = 0;
	// m,mp=12 t,u,v=201
	expCoeffs[292] = 0;
	// m,mp=12 t,u,v=011
	expCoeffs[293] = (-3*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=111
	expCoeffs[294] = 0;
	// m,mp=12 t,u,v=021
	expCoeffs[295] = -3/(8.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=002
	expCoeffs[296] = (3*(xPtoA + 2*xPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=102
	expCoeffs[297] = 0;
	// m,mp=12 t,u,v=012
	expCoeffs[298] = 0;
	// m,mp=12 t,u,v=003
	expCoeffs[299] = 3/(8.*std::pow(expGamma,3));
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff2And0(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {50};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=-20 t,u,v=000
	expCoeffs[0] = 6*xPtoA*yPtoA;
	// m,mp=-20 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=-20 t,u,v=200
	expCoeffs[2] = 0;
	// m,mp=-20 t,u,v=010
	expCoeffs[3] = (3*xPtoA)/expGamma;
	// m,mp=-20 t,u,v=110
	expCoeffs[4] = 0;
	// m,mp=-20 t,u,v=020
	expCoeffs[5] = 0;
	// m,mp=-20 t,u,v=001
	expCoeffs[6] = (3*yPtoA)/expGamma;
	// m,mp=-20 t,u,v=101
	expCoeffs[7] = 0;
	// m,mp=-20 t,u,v=011
	expCoeffs[8] = 3/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=002
	expCoeffs[9] = 0;
	// m,mp=-10 t,u,v=000
	expCoeffs[10] = 3*yPtoA*zPtoA;
	// m,mp=-10 t,u,v=100
	expCoeffs[11] = (3*yPtoA)/(2.*expGamma);
	// m,mp=-10 t,u,v=200
	expCoeffs[12] = 0;
	// m,mp=-10 t,u,v=010
	expCoeffs[13] = (3*zPtoA)/(2.*expGamma);
	// m,mp=-10 t,u,v=110
	expCoeffs[14] = 3/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=020
	expCoeffs[15] = 0;
	// m,mp=-10 t,u,v=001
	expCoeffs[16] = 0;
	// m,mp=-10 t,u,v=101
	expCoeffs[17] = 0;
	// m,mp=-10 t,u,v=011
	expCoeffs[18] = 0;
	// m,mp=-10 t,u,v=002
	expCoeffs[19] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[20] = (-d2PtoA + 3*std::pow(zPtoA,2))/2.;
	// m,mp=00 t,u,v=100
	expCoeffs[21] = zPtoA/expGamma;
	// m,mp=00 t,u,v=200
	expCoeffs[22] = 1/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=010
	expCoeffs[23] = -0.5*yPtoA/expGamma;
	// m,mp=00 t,u,v=110
	expCoeffs[24] = 0;
	// m,mp=00 t,u,v=020
	expCoeffs[25] = -0.125*1/std::pow(expGamma,2);
	// m,mp=00 t,u,v=001
	expCoeffs[26] = -0.5*xPtoA/expGamma;
	// m,mp=00 t,u,v=101
	expCoeffs[27] = 0;
	// m,mp=00 t,u,v=011
	expCoeffs[28] = 0;
	// m,mp=00 t,u,v=002
	expCoeffs[29] = -0.125*1/std::pow(expGamma,2);
	// m,mp=10 t,u,v=000
	expCoeffs[30] = 3*xPtoA*zPtoA;
	// m,mp=10 t,u,v=100
	expCoeffs[31] = (3*xPtoA)/(2.*expGamma);
	// m,mp=10 t,u,v=200
	expCoeffs[32] = 0;
	// m,mp=10 t,u,v=010
	expCoeffs[33] = 0;
	// m,mp=10 t,u,v=110
	expCoeffs[34] = 0;
	// m,mp=10 t,u,v=020
	expCoeffs[35] = 0;
	// m,mp=10 t,u,v=001
	expCoeffs[36] = (3*zPtoA)/(2.*expGamma);
	// m,mp=10 t,u,v=101
	expCoeffs[37] = 3/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=011
	expCoeffs[38] = 0;
	// m,mp=10 t,u,v=002
	expCoeffs[39] = 0;
	// m,mp=20 t,u,v=000
	expCoeffs[40] = 3*(xPtoA - yPtoA)*(xPtoA + yPtoA);
	// m,mp=20 t,u,v=100
	expCoeffs[41] = 0;
	// m,mp=20 t,u,v=200
	expCoeffs[42] = 0;
	// m,mp=20 t,u,v=010
	expCoeffs[43] = (-3*yPtoA)/expGamma;
	// m,mp=20 t,u,v=110
	expCoeffs[44] = 0;
	// m,mp=20 t,u,v=020
	expCoeffs[45] = -3/(4.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=001
	expCoeffs[46] = (3*xPtoA)/expGamma;
	// m,mp=20 t,u,v=101
	expCoeffs[47] = 0;
	// m,mp=20 t,u,v=011
	expCoeffs[48] = 0;
	// m,mp=20 t,u,v=002
	expCoeffs[49] = 3/(4.*std::pow(expGamma,2));
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff2And1(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {300};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=-2-1 t,u,v=000
	expCoeffs[0] = 3*xPtoA*(1/expGamma + 2*yPtoA*yPtoB);
	// m,mp=-2-1 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=-2-1 t,u,v=200
	expCoeffs[2] = 0;
	// m,mp=-2-1 t,u,v=300
	expCoeffs[3] = 0;
	// m,mp=-2-1 t,u,v=010
	expCoeffs[4] = (3*xPtoA*(yPtoA + yPtoB))/expGamma;
	// m,mp=-2-1 t,u,v=110
	expCoeffs[5] = 0;
	// m,mp=-2-1 t,u,v=210
	expCoeffs[6] = 0;
	// m,mp=-2-1 t,u,v=020
	expCoeffs[7] = (3*xPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=120
	expCoeffs[8] = 0;
	// m,mp=-2-1 t,u,v=030
	expCoeffs[9] = 0;
	// m,mp=-2-1 t,u,v=001
	expCoeffs[10] = (3 + 6*expGamma*yPtoA*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=101
	expCoeffs[11] = 0;
	// m,mp=-2-1 t,u,v=201
	expCoeffs[12] = 0;
	// m,mp=-2-1 t,u,v=011
	expCoeffs[13] = (3*(yPtoA + yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=111
	expCoeffs[14] = 0;
	// m,mp=-2-1 t,u,v=021
	expCoeffs[15] = 3/(4.*std::pow(expGamma,3));
	// m,mp=-2-1 t,u,v=002
	expCoeffs[16] = 0;
	// m,mp=-2-1 t,u,v=102
	expCoeffs[17] = 0;
	// m,mp=-2-1 t,u,v=012
	expCoeffs[18] = 0;
	// m,mp=-2-1 t,u,v=003
	expCoeffs[19] = 0;
	// m,mp=-20 t,u,v=000
	expCoeffs[20] = 6*xPtoA*yPtoA*zPtoB;
	// m,mp=-20 t,u,v=100
	expCoeffs[21] = (3*xPtoA*yPtoA)/expGamma;
	// m,mp=-20 t,u,v=200
	expCoeffs[22] = 0;
	// m,mp=-20 t,u,v=300
	expCoeffs[23] = 0;
	// m,mp=-20 t,u,v=010
	expCoeffs[24] = (3*xPtoA*zPtoB)/expGamma;
	// m,mp=-20 t,u,v=110
	expCoeffs[25] = (3*xPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=210
	expCoeffs[26] = 0;
	// m,mp=-20 t,u,v=020
	expCoeffs[27] = 0;
	// m,mp=-20 t,u,v=120
	expCoeffs[28] = 0;
	// m,mp=-20 t,u,v=030
	expCoeffs[29] = 0;
	// m,mp=-20 t,u,v=001
	expCoeffs[30] = (3*yPtoA*zPtoB)/expGamma;
	// m,mp=-20 t,u,v=101
	expCoeffs[31] = (3*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=201
	expCoeffs[32] = 0;
	// m,mp=-20 t,u,v=011
	expCoeffs[33] = (3*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=111
	expCoeffs[34] = 3/(4.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=021
	expCoeffs[35] = 0;
	// m,mp=-20 t,u,v=002
	expCoeffs[36] = 0;
	// m,mp=-20 t,u,v=102
	expCoeffs[37] = 0;
	// m,mp=-20 t,u,v=012
	expCoeffs[38] = 0;
	// m,mp=-20 t,u,v=003
	expCoeffs[39] = 0;
	// m,mp=-21 t,u,v=000
	expCoeffs[40] = 3*(1/expGamma + 2*xPtoA*xPtoB)*yPtoA;
	// m,mp=-21 t,u,v=100
	expCoeffs[41] = 0;
	// m,mp=-21 t,u,v=200
	expCoeffs[42] = 0;
	// m,mp=-21 t,u,v=300
	expCoeffs[43] = 0;
	// m,mp=-21 t,u,v=010
	expCoeffs[44] = (3 + 6*expGamma*xPtoA*xPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=110
	expCoeffs[45] = 0;
	// m,mp=-21 t,u,v=210
	expCoeffs[46] = 0;
	// m,mp=-21 t,u,v=020
	expCoeffs[47] = 0;
	// m,mp=-21 t,u,v=120
	expCoeffs[48] = 0;
	// m,mp=-21 t,u,v=030
	expCoeffs[49] = 0;
	// m,mp=-21 t,u,v=001
	expCoeffs[50] = (3*(xPtoA + xPtoB)*yPtoA)/expGamma;
	// m,mp=-21 t,u,v=101
	expCoeffs[51] = 0;
	// m,mp=-21 t,u,v=201
	expCoeffs[52] = 0;
	// m,mp=-21 t,u,v=011
	expCoeffs[53] = (3*(xPtoA + xPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=111
	expCoeffs[54] = 0;
	// m,mp=-21 t,u,v=021
	expCoeffs[55] = 0;
	// m,mp=-21 t,u,v=002
	expCoeffs[56] = (3*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=102
	expCoeffs[57] = 0;
	// m,mp=-21 t,u,v=012
	expCoeffs[58] = 3/(4.*std::pow(expGamma,3));
	// m,mp=-21 t,u,v=003
	expCoeffs[59] = 0;
	// m,mp=-1-1 t,u,v=000
	expCoeffs[60] = (3*(1/expGamma + 2*yPtoA*yPtoB)*zPtoA)/2.;
	// m,mp=-1-1 t,u,v=100
	expCoeffs[61] = (3 + 6*expGamma*yPtoA*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=200
	expCoeffs[62] = 0;
	// m,mp=-1-1 t,u,v=300
	expCoeffs[63] = 0;
	// m,mp=-1-1 t,u,v=010
	expCoeffs[64] = (3*(yPtoA + yPtoB)*zPtoA)/(2.*expGamma);
	// m,mp=-1-1 t,u,v=110
	expCoeffs[65] = (3*(yPtoA + yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=210
	expCoeffs[66] = 0;
	// m,mp=-1-1 t,u,v=020
	expCoeffs[67] = (3*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=120
	expCoeffs[68] = 3/(8.*std::pow(expGamma,3));
	// m,mp=-1-1 t,u,v=030
	expCoeffs[69] = 0;
	// m,mp=-1-1 t,u,v=001
	expCoeffs[70] = 0;
	// m,mp=-1-1 t,u,v=101
	expCoeffs[71] = 0;
	// m,mp=-1-1 t,u,v=201
	expCoeffs[72] = 0;
	// m,mp=-1-1 t,u,v=011
	expCoeffs[73] = 0;
	// m,mp=-1-1 t,u,v=111
	expCoeffs[74] = 0;
	// m,mp=-1-1 t,u,v=021
	expCoeffs[75] = 0;
	// m,mp=-1-1 t,u,v=002
	expCoeffs[76] = 0;
	// m,mp=-1-1 t,u,v=102
	expCoeffs[77] = 0;
	// m,mp=-1-1 t,u,v=012
	expCoeffs[78] = 0;
	// m,mp=-1-1 t,u,v=003
	expCoeffs[79] = 0;
	// m,mp=-10 t,u,v=000
	expCoeffs[80] = (3*yPtoA*(1/expGamma + 2*zPtoA*zPtoB))/2.;
	// m,mp=-10 t,u,v=100
	expCoeffs[81] = (3*yPtoA*(zPtoA + zPtoB))/(2.*expGamma);
	// m,mp=-10 t,u,v=200
	expCoeffs[82] = (3*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=300
	expCoeffs[83] = 0;
	// m,mp=-10 t,u,v=010
	expCoeffs[84] = (3 + 6*expGamma*zPtoA*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=110
	expCoeffs[85] = (3*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=210
	expCoeffs[86] = 3/(8.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=020
	expCoeffs[87] = 0;
	// m,mp=-10 t,u,v=120
	expCoeffs[88] = 0;
	// m,mp=-10 t,u,v=030
	expCoeffs[89] = 0;
	// m,mp=-10 t,u,v=001
	expCoeffs[90] = 0;
	// m,mp=-10 t,u,v=101
	expCoeffs[91] = 0;
	// m,mp=-10 t,u,v=201
	expCoeffs[92] = 0;
	// m,mp=-10 t,u,v=011
	expCoeffs[93] = 0;
	// m,mp=-10 t,u,v=111
	expCoeffs[94] = 0;
	// m,mp=-10 t,u,v=021
	expCoeffs[95] = 0;
	// m,mp=-10 t,u,v=002
	expCoeffs[96] = 0;
	// m,mp=-10 t,u,v=102
	expCoeffs[97] = 0;
	// m,mp=-10 t,u,v=012
	expCoeffs[98] = 0;
	// m,mp=-10 t,u,v=003
	expCoeffs[99] = 0;
	// m,mp=-11 t,u,v=000
	expCoeffs[100] = 3*xPtoB*yPtoA*zPtoA;
	// m,mp=-11 t,u,v=100
	expCoeffs[101] = (3*xPtoB*yPtoA)/(2.*expGamma);
	// m,mp=-11 t,u,v=200
	expCoeffs[102] = 0;
	// m,mp=-11 t,u,v=300
	expCoeffs[103] = 0;
	// m,mp=-11 t,u,v=010
	expCoeffs[104] = (3*xPtoB*zPtoA)/(2.*expGamma);
	// m,mp=-11 t,u,v=110
	expCoeffs[105] = (3*xPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=210
	expCoeffs[106] = 0;
	// m,mp=-11 t,u,v=020
	expCoeffs[107] = 0;
	// m,mp=-11 t,u,v=120
	expCoeffs[108] = 0;
	// m,mp=-11 t,u,v=030
	expCoeffs[109] = 0;
	// m,mp=-11 t,u,v=001
	expCoeffs[110] = (3*yPtoA*zPtoA)/(2.*expGamma);
	// m,mp=-11 t,u,v=101
	expCoeffs[111] = (3*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=201
	expCoeffs[112] = 0;
	// m,mp=-11 t,u,v=011
	expCoeffs[113] = (3*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=111
	expCoeffs[114] = 3/(8.*std::pow(expGamma,3));
	// m,mp=-11 t,u,v=021
	expCoeffs[115] = 0;
	// m,mp=-11 t,u,v=002
	expCoeffs[116] = 0;
	// m,mp=-11 t,u,v=102
	expCoeffs[117] = 0;
	// m,mp=-11 t,u,v=012
	expCoeffs[118] = 0;
	// m,mp=-11 t,u,v=003
	expCoeffs[119] = 0;
	// m,mp=0-1 t,u,v=000
	expCoeffs[120] = -0.5*(yPtoA + expGamma*yPtoB*(d2PtoA - 3*std::pow(zPtoA,2)))/expGamma;
	// m,mp=0-1 t,u,v=100
	expCoeffs[121] = (yPtoB*zPtoA)/expGamma;
	// m,mp=0-1 t,u,v=200
	expCoeffs[122] = yPtoB/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=300
	expCoeffs[123] = 0;
	// m,mp=0-1 t,u,v=010
	expCoeffs[124] = -0.25*(1 + expGamma*(d2PtoA + 2*yPtoA*yPtoB - 3*std::pow(zPtoA,2)))/std::pow(expGamma,2);
	// m,mp=0-1 t,u,v=110
	expCoeffs[125] = zPtoA/(2.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=210
	expCoeffs[126] = 1/(8.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=020
	expCoeffs[127] = -0.125*(2*yPtoA + yPtoB)/std::pow(expGamma,2);
	// m,mp=0-1 t,u,v=120
	expCoeffs[128] = 0;
	// m,mp=0-1 t,u,v=030
	expCoeffs[129] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=0-1 t,u,v=001
	expCoeffs[130] = -0.5*(xPtoA*yPtoB)/expGamma;
	// m,mp=0-1 t,u,v=101
	expCoeffs[131] = 0;
	// m,mp=0-1 t,u,v=201
	expCoeffs[132] = 0;
	// m,mp=0-1 t,u,v=011
	expCoeffs[133] = -0.25*xPtoA/std::pow(expGamma,2);
	// m,mp=0-1 t,u,v=111
	expCoeffs[134] = 0;
	// m,mp=0-1 t,u,v=021
	expCoeffs[135] = 0;
	// m,mp=0-1 t,u,v=002
	expCoeffs[136] = -0.125*yPtoB/std::pow(expGamma,2);
	// m,mp=0-1 t,u,v=102
	expCoeffs[137] = 0;
	// m,mp=0-1 t,u,v=012
	expCoeffs[138] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=0-1 t,u,v=003
	expCoeffs[139] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[140] = zPtoA/expGamma - ((d2PtoA - 3*std::pow(zPtoA,2))*zPtoB)/2.;
	// m,mp=00 t,u,v=100
	expCoeffs[141] = (2 - d2PtoA*expGamma + expGamma*zPtoA*(3*zPtoA + 4*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=200
	expCoeffs[142] = (2*zPtoA + zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=300
	expCoeffs[143] = 1/(8.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=010
	expCoeffs[144] = -0.5*(yPtoA*zPtoB)/expGamma;
	// m,mp=00 t,u,v=110
	expCoeffs[145] = -0.25*yPtoA/std::pow(expGamma,2);
	// m,mp=00 t,u,v=210
	expCoeffs[146] = 0;
	// m,mp=00 t,u,v=020
	expCoeffs[147] = -0.125*zPtoB/std::pow(expGamma,2);
	// m,mp=00 t,u,v=120
	expCoeffs[148] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=00 t,u,v=030
	expCoeffs[149] = 0;
	// m,mp=00 t,u,v=001
	expCoeffs[150] = -0.5*(xPtoA*zPtoB)/expGamma;
	// m,mp=00 t,u,v=101
	expCoeffs[151] = -0.25*xPtoA/std::pow(expGamma,2);
	// m,mp=00 t,u,v=201
	expCoeffs[152] = 0;
	// m,mp=00 t,u,v=011
	expCoeffs[153] = 0;
	// m,mp=00 t,u,v=111
	expCoeffs[154] = 0;
	// m,mp=00 t,u,v=021
	expCoeffs[155] = 0;
	// m,mp=00 t,u,v=002
	expCoeffs[156] = -0.125*zPtoB/std::pow(expGamma,2);
	// m,mp=00 t,u,v=102
	expCoeffs[157] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=00 t,u,v=012
	expCoeffs[158] = 0;
	// m,mp=00 t,u,v=003
	expCoeffs[159] = 0;
	// m,mp=01 t,u,v=000
	expCoeffs[160] = -0.5*(xPtoA + expGamma*xPtoB*(d2PtoA - 3*std::pow(zPtoA,2)))/expGamma;
	// m,mp=01 t,u,v=100
	expCoeffs[161] = (xPtoB*zPtoA)/expGamma;
	// m,mp=01 t,u,v=200
	expCoeffs[162] = xPtoB/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=300
	expCoeffs[163] = 0;
	// m,mp=01 t,u,v=010
	expCoeffs[164] = -0.5*(xPtoB*yPtoA)/expGamma;
	// m,mp=01 t,u,v=110
	expCoeffs[165] = 0;
	// m,mp=01 t,u,v=210
	expCoeffs[166] = 0;
	// m,mp=01 t,u,v=020
	expCoeffs[167] = -0.125*xPtoB/std::pow(expGamma,2);
	// m,mp=01 t,u,v=120
	expCoeffs[168] = 0;
	// m,mp=01 t,u,v=030
	expCoeffs[169] = 0;
	// m,mp=01 t,u,v=001
	expCoeffs[170] = -0.25*(1 + expGamma*(d2PtoA + 2*xPtoA*xPtoB - 3*std::pow(zPtoA,2)))/std::pow(expGamma,2);
	// m,mp=01 t,u,v=101
	expCoeffs[171] = zPtoA/(2.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=201
	expCoeffs[172] = 1/(8.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=011
	expCoeffs[173] = -0.25*yPtoA/std::pow(expGamma,2);
	// m,mp=01 t,u,v=111
	expCoeffs[174] = 0;
	// m,mp=01 t,u,v=021
	expCoeffs[175] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=01 t,u,v=002
	expCoeffs[176] = -0.125*(2*xPtoA + xPtoB)/std::pow(expGamma,2);
	// m,mp=01 t,u,v=102
	expCoeffs[177] = 0;
	// m,mp=01 t,u,v=012
	expCoeffs[178] = 0;
	// m,mp=01 t,u,v=003
	expCoeffs[179] = -0.0625*1/std::pow(expGamma,3);
	// m,mp=1-1 t,u,v=000
	expCoeffs[180] = 3*xPtoA*yPtoB*zPtoA;
	// m,mp=1-1 t,u,v=100
	expCoeffs[181] = (3*xPtoA*yPtoB)/(2.*expGamma);
	// m,mp=1-1 t,u,v=200
	expCoeffs[182] = 0;
	// m,mp=1-1 t,u,v=300
	expCoeffs[183] = 0;
	// m,mp=1-1 t,u,v=010
	expCoeffs[184] = (3*xPtoA*zPtoA)/(2.*expGamma);
	// m,mp=1-1 t,u,v=110
	expCoeffs[185] = (3*xPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=210
	expCoeffs[186] = 0;
	// m,mp=1-1 t,u,v=020
	expCoeffs[187] = 0;
	// m,mp=1-1 t,u,v=120
	expCoeffs[188] = 0;
	// m,mp=1-1 t,u,v=030
	expCoeffs[189] = 0;
	// m,mp=1-1 t,u,v=001
	expCoeffs[190] = (3*yPtoB*zPtoA)/(2.*expGamma);
	// m,mp=1-1 t,u,v=101
	expCoeffs[191] = (3*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=201
	expCoeffs[192] = 0;
	// m,mp=1-1 t,u,v=011
	expCoeffs[193] = (3*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=111
	expCoeffs[194] = 3/(8.*std::pow(expGamma,3));
	// m,mp=1-1 t,u,v=021
	expCoeffs[195] = 0;
	// m,mp=1-1 t,u,v=002
	expCoeffs[196] = 0;
	// m,mp=1-1 t,u,v=102
	expCoeffs[197] = 0;
	// m,mp=1-1 t,u,v=012
	expCoeffs[198] = 0;
	// m,mp=1-1 t,u,v=003
	expCoeffs[199] = 0;
	// m,mp=10 t,u,v=000
	expCoeffs[200] = (3*xPtoA*(1/expGamma + 2*zPtoA*zPtoB))/2.;
	// m,mp=10 t,u,v=100
	expCoeffs[201] = (3*xPtoA*(zPtoA + zPtoB))/(2.*expGamma);
	// m,mp=10 t,u,v=200
	expCoeffs[202] = (3*xPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=300
	expCoeffs[203] = 0;
	// m,mp=10 t,u,v=010
	expCoeffs[204] = 0;
	// m,mp=10 t,u,v=110
	expCoeffs[205] = 0;
	// m,mp=10 t,u,v=210
	expCoeffs[206] = 0;
	// m,mp=10 t,u,v=020
	expCoeffs[207] = 0;
	// m,mp=10 t,u,v=120
	expCoeffs[208] = 0;
	// m,mp=10 t,u,v=030
	expCoeffs[209] = 0;
	// m,mp=10 t,u,v=001
	expCoeffs[210] = (3 + 6*expGamma*zPtoA*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=101
	expCoeffs[211] = (3*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=201
	expCoeffs[212] = 3/(8.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=011
	expCoeffs[213] = 0;
	// m,mp=10 t,u,v=111
	expCoeffs[214] = 0;
	// m,mp=10 t,u,v=021
	expCoeffs[215] = 0;
	// m,mp=10 t,u,v=002
	expCoeffs[216] = 0;
	// m,mp=10 t,u,v=102
	expCoeffs[217] = 0;
	// m,mp=10 t,u,v=012
	expCoeffs[218] = 0;
	// m,mp=10 t,u,v=003
	expCoeffs[219] = 0;
	// m,mp=11 t,u,v=000
	expCoeffs[220] = (3*(1/expGamma + 2*xPtoA*xPtoB)*zPtoA)/2.;
	// m,mp=11 t,u,v=100
	expCoeffs[221] = (3 + 6*expGamma*xPtoA*xPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=200
	expCoeffs[222] = 0;
	// m,mp=11 t,u,v=300
	expCoeffs[223] = 0;
	// m,mp=11 t,u,v=010
	expCoeffs[224] = 0;
	// m,mp=11 t,u,v=110
	expCoeffs[225] = 0;
	// m,mp=11 t,u,v=210
	expCoeffs[226] = 0;
	// m,mp=11 t,u,v=020
	expCoeffs[227] = 0;
	// m,mp=11 t,u,v=120
	expCoeffs[228] = 0;
	// m,mp=11 t,u,v=030
	expCoeffs[229] = 0;
	// m,mp=11 t,u,v=001
	expCoeffs[230] = (3*(xPtoA + xPtoB)*zPtoA)/(2.*expGamma);
	// m,mp=11 t,u,v=101
	expCoeffs[231] = (3*(xPtoA + xPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=201
	expCoeffs[232] = 0;
	// m,mp=11 t,u,v=011
	expCoeffs[233] = 0;
	// m,mp=11 t,u,v=111
	expCoeffs[234] = 0;
	// m,mp=11 t,u,v=021
	expCoeffs[235] = 0;
	// m,mp=11 t,u,v=002
	expCoeffs[236] = (3*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=102
	expCoeffs[237] = 3/(8.*std::pow(expGamma,3));
	// m,mp=11 t,u,v=012
	expCoeffs[238] = 0;
	// m,mp=11 t,u,v=003
	expCoeffs[239] = 0;
	// m,mp=2-1 t,u,v=000
	expCoeffs[240] = 3*std::pow(xPtoA,2)*yPtoB - (3*yPtoA*(1 + expGamma*yPtoA*yPtoB))/expGamma;
	// m,mp=2-1 t,u,v=100
	expCoeffs[241] = 0;
	// m,mp=2-1 t,u,v=200
	expCoeffs[242] = 0;
	// m,mp=2-1 t,u,v=300
	expCoeffs[243] = 0;
	// m,mp=2-1 t,u,v=010
	expCoeffs[244] = (-3 + 3*expGamma*std::pow(xPtoA,2) - 3*expGamma*yPtoA*(yPtoA + 2*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=110
	expCoeffs[245] = 0;
	// m,mp=2-1 t,u,v=210
	expCoeffs[246] = 0;
	// m,mp=2-1 t,u,v=020
	expCoeffs[247] = (-3*(2*yPtoA + yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=120
	expCoeffs[248] = 0;
	// m,mp=2-1 t,u,v=030
	expCoeffs[249] = -3/(8.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=001
	expCoeffs[250] = (3*xPtoA*yPtoB)/expGamma;
	// m,mp=2-1 t,u,v=101
	expCoeffs[251] = 0;
	// m,mp=2-1 t,u,v=201
	expCoeffs[252] = 0;
	// m,mp=2-1 t,u,v=011
	expCoeffs[253] = (3*xPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=111
	expCoeffs[254] = 0;
	// m,mp=2-1 t,u,v=021
	expCoeffs[255] = 0;
	// m,mp=2-1 t,u,v=002
	expCoeffs[256] = (3*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=102
	expCoeffs[257] = 0;
	// m,mp=2-1 t,u,v=012
	expCoeffs[258] = 3/(8.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=003
	expCoeffs[259] = 0;
	// m,mp=20 t,u,v=000
	expCoeffs[260] = 3*(xPtoA - yPtoA)*(xPtoA + yPtoA)*zPtoB;
	// m,mp=20 t,u,v=100
	expCoeffs[261] = (3*(xPtoA - yPtoA)*(xPtoA + yPtoA))/(2.*expGamma);
	// m,mp=20 t,u,v=200
	expCoeffs[262] = 0;
	// m,mp=20 t,u,v=300
	expCoeffs[263] = 0;
	// m,mp=20 t,u,v=010
	expCoeffs[264] = (-3*yPtoA*zPtoB)/expGamma;
	// m,mp=20 t,u,v=110
	expCoeffs[265] = (-3*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=210
	expCoeffs[266] = 0;
	// m,mp=20 t,u,v=020
	expCoeffs[267] = (-3*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=120
	expCoeffs[268] = -3/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=030
	expCoeffs[269] = 0;
	// m,mp=20 t,u,v=001
	expCoeffs[270] = (3*xPtoA*zPtoB)/expGamma;
	// m,mp=20 t,u,v=101
	expCoeffs[271] = (3*xPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=201
	expCoeffs[272] = 0;
	// m,mp=20 t,u,v=011
	expCoeffs[273] = 0;
	// m,mp=20 t,u,v=111
	expCoeffs[274] = 0;
	// m,mp=20 t,u,v=021
	expCoeffs[275] = 0;
	// m,mp=20 t,u,v=002
	expCoeffs[276] = (3*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=102
	expCoeffs[277] = 3/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=012
	expCoeffs[278] = 0;
	// m,mp=20 t,u,v=003
	expCoeffs[279] = 0;
	// m,mp=21 t,u,v=000
	expCoeffs[280] = (3*(xPtoA + expGamma*xPtoB*(xPtoA - yPtoA)*(xPtoA + yPtoA)))/expGamma;
	// m,mp=21 t,u,v=100
	expCoeffs[281] = 0;
	// m,mp=21 t,u,v=200
	expCoeffs[282] = 0;
	// m,mp=21 t,u,v=300
	expCoeffs[283] = 0;
	// m,mp=21 t,u,v=010
	expCoeffs[284] = (-3*xPtoB*yPtoA)/expGamma;
	// m,mp=21 t,u,v=110
	expCoeffs[285] = 0;
	// m,mp=21 t,u,v=210
	expCoeffs[286] = 0;
	// m,mp=21 t,u,v=020
	expCoeffs[287] = (-3*xPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=120
	expCoeffs[288] = 0;
	// m,mp=21 t,u,v=030
	expCoeffs[289] = 0;
	// m,mp=21 t,u,v=001
	expCoeffs[290] = (3*(1 + expGamma*(std::pow(xPtoA,2) + 2*xPtoA*xPtoB - std::pow(yPtoA,2))))/(2.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=101
	expCoeffs[291] = 0;
	// m,mp=21 t,u,v=201
	expCoeffs[292] = 0;
	// m,mp=21 t,u,v=011
	expCoeffs[293] = (-3*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=111
	expCoeffs[294] = 0;
	// m,mp=21 t,u,v=021
	expCoeffs[295] = -3/(8.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=002
	expCoeffs[296] = (3*(2*xPtoA + xPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=102
	expCoeffs[297] = 0;
	// m,mp=21 t,u,v=012
	expCoeffs[298] = 0;
	// m,mp=21 t,u,v=003
	expCoeffs[299] = 3/(8.*std::pow(expGamma,3));
	return expCoeffs;
}

#pragma acc routine seq
std::vector<double> ECoeff2And2(double expAlpha, double expBeta,
	std::array<double,3> muCoords, std::array<double,3> nuCoords)
{
	double expGamma{expAlpha + expBeta};
	double betaOverGamma{expBeta / expGamma};
	double xPtoA {betaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoA {betaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoA {betaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoA {xPtoA*xPtoA + yPtoA*yPtoA + zPtoA*zPtoA};
	double alphaOverGamma{-expAlpha / expGamma};
	double xPtoB {alphaOverGamma*(nuCoords[0] - muCoords[0])};
	double yPtoB {alphaOverGamma*(nuCoords[1] - muCoords[1])};
	double zPtoB {alphaOverGamma*(nuCoords[2] - muCoords[2])};
	double d2PtoB {xPtoB*xPtoB + yPtoB*yPtoB + zPtoB*zPtoB};
	int nCoeffs {875};
	std::vector<double> expCoeffs(nCoeffs,0);
	// m,mp=-2-2 t,u,v=000
	expCoeffs[0] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*(1 + 2*expGamma*yPtoA*yPtoB))/std::pow(expGamma,2);
	// m,mp=-2-2 t,u,v=100
	expCoeffs[1] = 0;
	// m,mp=-2-2 t,u,v=200
	expCoeffs[2] = 0;
	// m,mp=-2-2 t,u,v=300
	expCoeffs[3] = 0;
	// m,mp=-2-2 t,u,v=400
	expCoeffs[4] = 0;
	// m,mp=-2-2 t,u,v=010
	expCoeffs[5] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*(yPtoA + yPtoB))/std::pow(expGamma,2);
	// m,mp=-2-2 t,u,v=110
	expCoeffs[6] = 0;
	// m,mp=-2-2 t,u,v=210
	expCoeffs[7] = 0;
	// m,mp=-2-2 t,u,v=310
	expCoeffs[8] = 0;
	// m,mp=-2-2 t,u,v=020
	expCoeffs[9] = (9*(1 + 2*expGamma*xPtoA*xPtoB))/(2.*std::pow(expGamma,3));
	// m,mp=-2-2 t,u,v=120
	expCoeffs[10] = 0;
	// m,mp=-2-2 t,u,v=220
	expCoeffs[11] = 0;
	// m,mp=-2-2 t,u,v=030
	expCoeffs[12] = 0;
	// m,mp=-2-2 t,u,v=130
	expCoeffs[13] = 0;
	// m,mp=-2-2 t,u,v=040
	expCoeffs[14] = 0;
	// m,mp=-2-2 t,u,v=001
	expCoeffs[15] = (9*(xPtoA + xPtoB)*(1 + 2*expGamma*yPtoA*yPtoB))/std::pow(expGamma,2);
	// m,mp=-2-2 t,u,v=101
	expCoeffs[16] = 0;
	// m,mp=-2-2 t,u,v=201
	expCoeffs[17] = 0;
	// m,mp=-2-2 t,u,v=301
	expCoeffs[18] = 0;
	// m,mp=-2-2 t,u,v=011
	expCoeffs[19] = (9*(xPtoA + xPtoB)*(yPtoA + yPtoB))/std::pow(expGamma,2);
	// m,mp=-2-2 t,u,v=111
	expCoeffs[20] = 0;
	// m,mp=-2-2 t,u,v=211
	expCoeffs[21] = 0;
	// m,mp=-2-2 t,u,v=021
	expCoeffs[22] = (9*(xPtoA + xPtoB))/(2.*std::pow(expGamma,3));
	// m,mp=-2-2 t,u,v=121
	expCoeffs[23] = 0;
	// m,mp=-2-2 t,u,v=031
	expCoeffs[24] = 0;
	// m,mp=-2-2 t,u,v=002
	expCoeffs[25] = (9*(1 + 2*expGamma*yPtoA*yPtoB))/(2.*std::pow(expGamma,3));
	// m,mp=-2-2 t,u,v=102
	expCoeffs[26] = 0;
	// m,mp=-2-2 t,u,v=202
	expCoeffs[27] = 0;
	// m,mp=-2-2 t,u,v=012
	expCoeffs[28] = (9*(yPtoA + yPtoB))/(2.*std::pow(expGamma,3));
	// m,mp=-2-2 t,u,v=112
	expCoeffs[29] = 0;
	// m,mp=-2-2 t,u,v=022
	expCoeffs[30] = 9/(4.*std::pow(expGamma,4));
	// m,mp=-2-2 t,u,v=003
	expCoeffs[31] = 0;
	// m,mp=-2-2 t,u,v=103
	expCoeffs[32] = 0;
	// m,mp=-2-2 t,u,v=013
	expCoeffs[33] = 0;
	// m,mp=-2-2 t,u,v=004
	expCoeffs[34] = 0;
	// m,mp=-2-1 t,u,v=000
	expCoeffs[35] = (9*xPtoA*(1 + 2*expGamma*yPtoA*yPtoB)*zPtoB)/expGamma;
	// m,mp=-2-1 t,u,v=100
	expCoeffs[36] = (9*(xPtoA + 2*expGamma*xPtoA*yPtoA*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=200
	expCoeffs[37] = 0;
	// m,mp=-2-1 t,u,v=300
	expCoeffs[38] = 0;
	// m,mp=-2-1 t,u,v=400
	expCoeffs[39] = 0;
	// m,mp=-2-1 t,u,v=010
	expCoeffs[40] = (9*xPtoA*(yPtoA + yPtoB)*zPtoB)/expGamma;
	// m,mp=-2-1 t,u,v=110
	expCoeffs[41] = (9*xPtoA*(yPtoA + yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=210
	expCoeffs[42] = 0;
	// m,mp=-2-1 t,u,v=310
	expCoeffs[43] = 0;
	// m,mp=-2-1 t,u,v=020
	expCoeffs[44] = (9*xPtoA*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=120
	expCoeffs[45] = (9*xPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-2-1 t,u,v=220
	expCoeffs[46] = 0;
	// m,mp=-2-1 t,u,v=030
	expCoeffs[47] = 0;
	// m,mp=-2-1 t,u,v=130
	expCoeffs[48] = 0;
	// m,mp=-2-1 t,u,v=040
	expCoeffs[49] = 0;
	// m,mp=-2-1 t,u,v=001
	expCoeffs[50] = (9*(1 + 2*expGamma*yPtoA*yPtoB)*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=101
	expCoeffs[51] = (9*(1 + 2*expGamma*yPtoA*yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-2-1 t,u,v=201
	expCoeffs[52] = 0;
	// m,mp=-2-1 t,u,v=301
	expCoeffs[53] = 0;
	// m,mp=-2-1 t,u,v=011
	expCoeffs[54] = (9*(yPtoA + yPtoB)*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-2-1 t,u,v=111
	expCoeffs[55] = (9*(yPtoA + yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-2-1 t,u,v=211
	expCoeffs[56] = 0;
	// m,mp=-2-1 t,u,v=021
	expCoeffs[57] = (9*zPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=-2-1 t,u,v=121
	expCoeffs[58] = 9/(8.*std::pow(expGamma,4));
	// m,mp=-2-1 t,u,v=031
	expCoeffs[59] = 0;
	// m,mp=-2-1 t,u,v=002
	expCoeffs[60] = 0;
	// m,mp=-2-1 t,u,v=102
	expCoeffs[61] = 0;
	// m,mp=-2-1 t,u,v=202
	expCoeffs[62] = 0;
	// m,mp=-2-1 t,u,v=012
	expCoeffs[63] = 0;
	// m,mp=-2-1 t,u,v=112
	expCoeffs[64] = 0;
	// m,mp=-2-1 t,u,v=022
	expCoeffs[65] = 0;
	// m,mp=-2-1 t,u,v=003
	expCoeffs[66] = 0;
	// m,mp=-2-1 t,u,v=103
	expCoeffs[67] = 0;
	// m,mp=-2-1 t,u,v=013
	expCoeffs[68] = 0;
	// m,mp=-2-1 t,u,v=004
	expCoeffs[69] = 0;
	// m,mp=-20 t,u,v=000
	expCoeffs[70] = (-3*(xPtoB*yPtoA + xPtoA*(yPtoB + expGamma*yPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))))/expGamma;
	// m,mp=-20 t,u,v=100
	expCoeffs[71] = (6*xPtoA*yPtoA*zPtoB)/expGamma;
	// m,mp=-20 t,u,v=200
	expCoeffs[72] = (3*xPtoA*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=300
	expCoeffs[73] = 0;
	// m,mp=-20 t,u,v=400
	expCoeffs[74] = 0;
	// m,mp=-20 t,u,v=010
	expCoeffs[75] = (-3*(xPtoA + xPtoB + expGamma*xPtoA*(d2PtoB + 2*yPtoA*yPtoB - 3*std::pow(zPtoB,2))))/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=110
	expCoeffs[76] = (3*xPtoA*zPtoB)/std::pow(expGamma,2);
	// m,mp=-20 t,u,v=210
	expCoeffs[77] = (3*xPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=310
	expCoeffs[78] = 0;
	// m,mp=-20 t,u,v=020
	expCoeffs[79] = (-3*xPtoA*(yPtoA + 2*yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=120
	expCoeffs[80] = 0;
	// m,mp=-20 t,u,v=220
	expCoeffs[81] = 0;
	// m,mp=-20 t,u,v=030
	expCoeffs[82] = (-3*xPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=130
	expCoeffs[83] = 0;
	// m,mp=-20 t,u,v=040
	expCoeffs[84] = 0;
	// m,mp=-20 t,u,v=001
	expCoeffs[85] = (-3*(yPtoA + yPtoB + expGamma*yPtoA*(d2PtoB + 2*xPtoA*xPtoB - 3*std::pow(zPtoB,2))))/(2.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=101
	expCoeffs[86] = (3*yPtoA*zPtoB)/std::pow(expGamma,2);
	// m,mp=-20 t,u,v=201
	expCoeffs[87] = (3*yPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=301
	expCoeffs[88] = 0;
	// m,mp=-20 t,u,v=011
	expCoeffs[89] = (-3*(2 + expGamma*(d2PtoB + 2*xPtoA*xPtoB + 2*yPtoA*yPtoB - 3*std::pow(zPtoB,2))))/(4.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=111
	expCoeffs[90] = (3*zPtoB)/(2.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=211
	expCoeffs[91] = 3/(8.*std::pow(expGamma,4));
	// m,mp=-20 t,u,v=021
	expCoeffs[92] = (-3*(yPtoA + 2*yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=121
	expCoeffs[93] = 0;
	// m,mp=-20 t,u,v=031
	expCoeffs[94] = -3/(16.*std::pow(expGamma,4));
	// m,mp=-20 t,u,v=002
	expCoeffs[95] = (-3*(xPtoA + 2*xPtoB)*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-20 t,u,v=102
	expCoeffs[96] = 0;
	// m,mp=-20 t,u,v=202
	expCoeffs[97] = 0;
	// m,mp=-20 t,u,v=012
	expCoeffs[98] = (-3*(xPtoA + 2*xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=112
	expCoeffs[99] = 0;
	// m,mp=-20 t,u,v=022
	expCoeffs[100] = 0;
	// m,mp=-20 t,u,v=003
	expCoeffs[101] = (-3*yPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-20 t,u,v=103
	expCoeffs[102] = 0;
	// m,mp=-20 t,u,v=013
	expCoeffs[103] = -3/(16.*std::pow(expGamma,4));
	// m,mp=-20 t,u,v=004
	expCoeffs[104] = 0;
	// m,mp=-21 t,u,v=000
	expCoeffs[105] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*yPtoA*zPtoB)/expGamma;
	// m,mp=-21 t,u,v=100
	expCoeffs[106] = (9*(yPtoA + 2*expGamma*xPtoA*xPtoB*yPtoA))/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=200
	expCoeffs[107] = 0;
	// m,mp=-21 t,u,v=300
	expCoeffs[108] = 0;
	// m,mp=-21 t,u,v=400
	expCoeffs[109] = 0;
	// m,mp=-21 t,u,v=010
	expCoeffs[110] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=110
	expCoeffs[111] = (9*(1 + 2*expGamma*xPtoA*xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-21 t,u,v=210
	expCoeffs[112] = 0;
	// m,mp=-21 t,u,v=310
	expCoeffs[113] = 0;
	// m,mp=-21 t,u,v=020
	expCoeffs[114] = 0;
	// m,mp=-21 t,u,v=120
	expCoeffs[115] = 0;
	// m,mp=-21 t,u,v=220
	expCoeffs[116] = 0;
	// m,mp=-21 t,u,v=030
	expCoeffs[117] = 0;
	// m,mp=-21 t,u,v=130
	expCoeffs[118] = 0;
	// m,mp=-21 t,u,v=040
	expCoeffs[119] = 0;
	// m,mp=-21 t,u,v=001
	expCoeffs[120] = (9*(xPtoA + xPtoB)*yPtoA*zPtoB)/expGamma;
	// m,mp=-21 t,u,v=101
	expCoeffs[121] = (9*(xPtoA + xPtoB)*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=201
	expCoeffs[122] = 0;
	// m,mp=-21 t,u,v=301
	expCoeffs[123] = 0;
	// m,mp=-21 t,u,v=011
	expCoeffs[124] = (9*(xPtoA + xPtoB)*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=111
	expCoeffs[125] = (9*(xPtoA + xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-21 t,u,v=211
	expCoeffs[126] = 0;
	// m,mp=-21 t,u,v=021
	expCoeffs[127] = 0;
	// m,mp=-21 t,u,v=121
	expCoeffs[128] = 0;
	// m,mp=-21 t,u,v=031
	expCoeffs[129] = 0;
	// m,mp=-21 t,u,v=002
	expCoeffs[130] = (9*yPtoA*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=-21 t,u,v=102
	expCoeffs[131] = (9*yPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-21 t,u,v=202
	expCoeffs[132] = 0;
	// m,mp=-21 t,u,v=012
	expCoeffs[133] = (9*zPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=-21 t,u,v=112
	expCoeffs[134] = 9/(8.*std::pow(expGamma,4));
	// m,mp=-21 t,u,v=022
	expCoeffs[135] = 0;
	// m,mp=-21 t,u,v=003
	expCoeffs[136] = 0;
	// m,mp=-21 t,u,v=103
	expCoeffs[137] = 0;
	// m,mp=-21 t,u,v=013
	expCoeffs[138] = 0;
	// m,mp=-21 t,u,v=004
	expCoeffs[139] = 0;
	// m,mp=-22 t,u,v=000
	expCoeffs[140] = (18*(xPtoB*yPtoA + expGamma*xPtoA*std::pow(xPtoB,2)*yPtoA - xPtoA*yPtoB*(1 + expGamma*yPtoA*yPtoB)))/expGamma;
	// m,mp=-22 t,u,v=100
	expCoeffs[141] = 0;
	// m,mp=-22 t,u,v=200
	expCoeffs[142] = 0;
	// m,mp=-22 t,u,v=300
	expCoeffs[143] = 0;
	// m,mp=-22 t,u,v=400
	expCoeffs[144] = 0;
	// m,mp=-22 t,u,v=010
	expCoeffs[145] = (9*(xPtoB + xPtoA*(-1 + expGamma*std::pow(xPtoB,2) - expGamma*yPtoB*(2*yPtoA + yPtoB))))/std::pow(expGamma,2);
	// m,mp=-22 t,u,v=110
	expCoeffs[146] = 0;
	// m,mp=-22 t,u,v=210
	expCoeffs[147] = 0;
	// m,mp=-22 t,u,v=310
	expCoeffs[148] = 0;
	// m,mp=-22 t,u,v=020
	expCoeffs[149] = (-9*xPtoA*(yPtoA + 2*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-22 t,u,v=120
	expCoeffs[150] = 0;
	// m,mp=-22 t,u,v=220
	expCoeffs[151] = 0;
	// m,mp=-22 t,u,v=030
	expCoeffs[152] = (-9*xPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-22 t,u,v=130
	expCoeffs[153] = 0;
	// m,mp=-22 t,u,v=040
	expCoeffs[154] = 0;
	// m,mp=-22 t,u,v=001
	expCoeffs[155] = (9*(yPtoA - yPtoB + expGamma*yPtoA*(2*xPtoA*xPtoB + std::pow(xPtoB,2) - std::pow(yPtoB,2))))/std::pow(expGamma,2);
	// m,mp=-22 t,u,v=101
	expCoeffs[156] = 0;
	// m,mp=-22 t,u,v=201
	expCoeffs[157] = 0;
	// m,mp=-22 t,u,v=301
	expCoeffs[158] = 0;
	// m,mp=-22 t,u,v=011
	expCoeffs[159] = (9*(2*xPtoA*xPtoB + std::pow(xPtoB,2) - yPtoB*(2*yPtoA + yPtoB)))/(2.*std::pow(expGamma,2));
	// m,mp=-22 t,u,v=111
	expCoeffs[160] = 0;
	// m,mp=-22 t,u,v=211
	expCoeffs[161] = 0;
	// m,mp=-22 t,u,v=021
	expCoeffs[162] = (-9*(yPtoA + 2*yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-22 t,u,v=121
	expCoeffs[163] = 0;
	// m,mp=-22 t,u,v=031
	expCoeffs[164] = -9/(8.*std::pow(expGamma,4));
	// m,mp=-22 t,u,v=002
	expCoeffs[165] = (9*(xPtoA + 2*xPtoB)*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-22 t,u,v=102
	expCoeffs[166] = 0;
	// m,mp=-22 t,u,v=202
	expCoeffs[167] = 0;
	// m,mp=-22 t,u,v=012
	expCoeffs[168] = (9*(xPtoA + 2*xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-22 t,u,v=112
	expCoeffs[169] = 0;
	// m,mp=-22 t,u,v=022
	expCoeffs[170] = 0;
	// m,mp=-22 t,u,v=003
	expCoeffs[171] = (9*yPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-22 t,u,v=103
	expCoeffs[172] = 0;
	// m,mp=-22 t,u,v=013
	expCoeffs[173] = 9/(8.*std::pow(expGamma,4));
	// m,mp=-22 t,u,v=004
	expCoeffs[174] = 0;
	// m,mp=-1-2 t,u,v=000
	expCoeffs[175] = 9*xPtoB*(1/expGamma + 2*yPtoA*yPtoB)*zPtoA;
	// m,mp=-1-2 t,u,v=100
	expCoeffs[176] = (9*xPtoB*(1 + 2*expGamma*yPtoA*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=200
	expCoeffs[177] = 0;
	// m,mp=-1-2 t,u,v=300
	expCoeffs[178] = 0;
	// m,mp=-1-2 t,u,v=400
	expCoeffs[179] = 0;
	// m,mp=-1-2 t,u,v=010
	expCoeffs[180] = (9*xPtoB*(yPtoA + yPtoB)*zPtoA)/expGamma;
	// m,mp=-1-2 t,u,v=110
	expCoeffs[181] = (9*xPtoB*(yPtoA + yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=210
	expCoeffs[182] = 0;
	// m,mp=-1-2 t,u,v=310
	expCoeffs[183] = 0;
	// m,mp=-1-2 t,u,v=020
	expCoeffs[184] = (9*xPtoB*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=120
	expCoeffs[185] = (9*xPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=-1-2 t,u,v=220
	expCoeffs[186] = 0;
	// m,mp=-1-2 t,u,v=030
	expCoeffs[187] = 0;
	// m,mp=-1-2 t,u,v=130
	expCoeffs[188] = 0;
	// m,mp=-1-2 t,u,v=040
	expCoeffs[189] = 0;
	// m,mp=-1-2 t,u,v=001
	expCoeffs[190] = (9*(1 + 2*expGamma*yPtoA*yPtoB)*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=101
	expCoeffs[191] = (9*(1 + 2*expGamma*yPtoA*yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-1-2 t,u,v=201
	expCoeffs[192] = 0;
	// m,mp=-1-2 t,u,v=301
	expCoeffs[193] = 0;
	// m,mp=-1-2 t,u,v=011
	expCoeffs[194] = (9*(yPtoA + yPtoB)*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-1-2 t,u,v=111
	expCoeffs[195] = (9*(yPtoA + yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=-1-2 t,u,v=211
	expCoeffs[196] = 0;
	// m,mp=-1-2 t,u,v=021
	expCoeffs[197] = (9*zPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=-1-2 t,u,v=121
	expCoeffs[198] = 9/(8.*std::pow(expGamma,4));
	// m,mp=-1-2 t,u,v=031
	expCoeffs[199] = 0;
	// m,mp=-1-2 t,u,v=002
	expCoeffs[200] = 0;
	// m,mp=-1-2 t,u,v=102
	expCoeffs[201] = 0;
	// m,mp=-1-2 t,u,v=202
	expCoeffs[202] = 0;
	// m,mp=-1-2 t,u,v=012
	expCoeffs[203] = 0;
	// m,mp=-1-2 t,u,v=112
	expCoeffs[204] = 0;
	// m,mp=-1-2 t,u,v=022
	expCoeffs[205] = 0;
	// m,mp=-1-2 t,u,v=003
	expCoeffs[206] = 0;
	// m,mp=-1-2 t,u,v=103
	expCoeffs[207] = 0;
	// m,mp=-1-2 t,u,v=013
	expCoeffs[208] = 0;
	// m,mp=-1-2 t,u,v=004
	expCoeffs[209] = 0;
	// m,mp=-1-1 t,u,v=000
	expCoeffs[210] = (9*(1 + 2*expGamma*yPtoA*yPtoB)*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=100
	expCoeffs[211] = (9*(1 + 2*expGamma*yPtoA*yPtoB)*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=200
	expCoeffs[212] = (9*(1 + 2*expGamma*yPtoA*yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-1-1 t,u,v=300
	expCoeffs[213] = 0;
	// m,mp=-1-1 t,u,v=400
	expCoeffs[214] = 0;
	// m,mp=-1-1 t,u,v=010
	expCoeffs[215] = (9*(yPtoA + yPtoB)*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=110
	expCoeffs[216] = (9*(yPtoA + yPtoB)*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-1-1 t,u,v=210
	expCoeffs[217] = (9*(yPtoA + yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-1-1 t,u,v=310
	expCoeffs[218] = 0;
	// m,mp=-1-1 t,u,v=020
	expCoeffs[219] = (9*(1 + 2*expGamma*zPtoA*zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-1-1 t,u,v=120
	expCoeffs[220] = (9*(zPtoA + zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-1-1 t,u,v=220
	expCoeffs[221] = 9/(16.*std::pow(expGamma,4));
	// m,mp=-1-1 t,u,v=030
	expCoeffs[222] = 0;
	// m,mp=-1-1 t,u,v=130
	expCoeffs[223] = 0;
	// m,mp=-1-1 t,u,v=040
	expCoeffs[224] = 0;
	// m,mp=-1-1 t,u,v=001
	expCoeffs[225] = 0;
	// m,mp=-1-1 t,u,v=101
	expCoeffs[226] = 0;
	// m,mp=-1-1 t,u,v=201
	expCoeffs[227] = 0;
	// m,mp=-1-1 t,u,v=301
	expCoeffs[228] = 0;
	// m,mp=-1-1 t,u,v=011
	expCoeffs[229] = 0;
	// m,mp=-1-1 t,u,v=111
	expCoeffs[230] = 0;
	// m,mp=-1-1 t,u,v=211
	expCoeffs[231] = 0;
	// m,mp=-1-1 t,u,v=021
	expCoeffs[232] = 0;
	// m,mp=-1-1 t,u,v=121
	expCoeffs[233] = 0;
	// m,mp=-1-1 t,u,v=031
	expCoeffs[234] = 0;
	// m,mp=-1-1 t,u,v=002
	expCoeffs[235] = 0;
	// m,mp=-1-1 t,u,v=102
	expCoeffs[236] = 0;
	// m,mp=-1-1 t,u,v=202
	expCoeffs[237] = 0;
	// m,mp=-1-1 t,u,v=012
	expCoeffs[238] = 0;
	// m,mp=-1-1 t,u,v=112
	expCoeffs[239] = 0;
	// m,mp=-1-1 t,u,v=022
	expCoeffs[240] = 0;
	// m,mp=-1-1 t,u,v=003
	expCoeffs[241] = 0;
	// m,mp=-1-1 t,u,v=103
	expCoeffs[242] = 0;
	// m,mp=-1-1 t,u,v=013
	expCoeffs[243] = 0;
	// m,mp=-1-1 t,u,v=004
	expCoeffs[244] = 0;
	// m,mp=-10 t,u,v=000
	expCoeffs[245] = (3*(2*yPtoA*zPtoB - zPtoA*(yPtoB + expGamma*yPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))))/(2.*expGamma);
	// m,mp=-10 t,u,v=100
	expCoeffs[246] = (-3*(yPtoB + yPtoA*(-2 + d2PtoB*expGamma - expGamma*zPtoB*(4*zPtoA + 3*zPtoB))))/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=200
	expCoeffs[247] = (3*yPtoA*(zPtoA + 2*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=300
	expCoeffs[248] = (3*yPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=400
	expCoeffs[249] = 0;
	// m,mp=-10 t,u,v=010
	expCoeffs[250] = (-3*(zPtoA - 2*zPtoB + expGamma*zPtoA*(d2PtoB + 2*yPtoA*yPtoB - 3*std::pow(zPtoB,2))))/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=110
	expCoeffs[251] = (3 - 3*expGamma*(d2PtoB + 2*yPtoA*yPtoB - zPtoB*(4*zPtoA + 3*zPtoB)))/(8.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=210
	expCoeffs[252] = (3*(zPtoA + 2*zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=310
	expCoeffs[253] = 3/(16.*std::pow(expGamma,4));
	// m,mp=-10 t,u,v=020
	expCoeffs[254] = (-3*(yPtoA + 2*yPtoB)*zPtoA)/(8.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=120
	expCoeffs[255] = (-3*(yPtoA + 2*yPtoB))/(16.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=220
	expCoeffs[256] = 0;
	// m,mp=-10 t,u,v=030
	expCoeffs[257] = (-3*zPtoA)/(16.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=130
	expCoeffs[258] = -3/(32.*std::pow(expGamma,4));
	// m,mp=-10 t,u,v=040
	expCoeffs[259] = 0;
	// m,mp=-10 t,u,v=001
	expCoeffs[260] = (-3*xPtoB*yPtoA*zPtoA)/(2.*expGamma);
	// m,mp=-10 t,u,v=101
	expCoeffs[261] = (-3*xPtoB*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=201
	expCoeffs[262] = 0;
	// m,mp=-10 t,u,v=301
	expCoeffs[263] = 0;
	// m,mp=-10 t,u,v=011
	expCoeffs[264] = (-3*xPtoB*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=111
	expCoeffs[265] = (-3*xPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=211
	expCoeffs[266] = 0;
	// m,mp=-10 t,u,v=021
	expCoeffs[267] = 0;
	// m,mp=-10 t,u,v=121
	expCoeffs[268] = 0;
	// m,mp=-10 t,u,v=031
	expCoeffs[269] = 0;
	// m,mp=-10 t,u,v=002
	expCoeffs[270] = (-3*yPtoA*zPtoA)/(8.*std::pow(expGamma,2));
	// m,mp=-10 t,u,v=102
	expCoeffs[271] = (-3*yPtoA)/(16.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=202
	expCoeffs[272] = 0;
	// m,mp=-10 t,u,v=012
	expCoeffs[273] = (-3*zPtoA)/(16.*std::pow(expGamma,3));
	// m,mp=-10 t,u,v=112
	expCoeffs[274] = -3/(32.*std::pow(expGamma,4));
	// m,mp=-10 t,u,v=022
	expCoeffs[275] = 0;
	// m,mp=-10 t,u,v=003
	expCoeffs[276] = 0;
	// m,mp=-10 t,u,v=103
	expCoeffs[277] = 0;
	// m,mp=-10 t,u,v=013
	expCoeffs[278] = 0;
	// m,mp=-10 t,u,v=004
	expCoeffs[279] = 0;
	// m,mp=-11 t,u,v=000
	expCoeffs[280] = (9*xPtoB*yPtoA*(1/expGamma + 2*zPtoA*zPtoB))/2.;
	// m,mp=-11 t,u,v=100
	expCoeffs[281] = (9*xPtoB*yPtoA*(zPtoA + zPtoB))/(2.*expGamma);
	// m,mp=-11 t,u,v=200
	expCoeffs[282] = (9*xPtoB*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=300
	expCoeffs[283] = 0;
	// m,mp=-11 t,u,v=400
	expCoeffs[284] = 0;
	// m,mp=-11 t,u,v=010
	expCoeffs[285] = (9*xPtoB*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=110
	expCoeffs[286] = (9*xPtoB*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=210
	expCoeffs[287] = (9*xPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=-11 t,u,v=310
	expCoeffs[288] = 0;
	// m,mp=-11 t,u,v=020
	expCoeffs[289] = 0;
	// m,mp=-11 t,u,v=120
	expCoeffs[290] = 0;
	// m,mp=-11 t,u,v=220
	expCoeffs[291] = 0;
	// m,mp=-11 t,u,v=030
	expCoeffs[292] = 0;
	// m,mp=-11 t,u,v=130
	expCoeffs[293] = 0;
	// m,mp=-11 t,u,v=040
	expCoeffs[294] = 0;
	// m,mp=-11 t,u,v=001
	expCoeffs[295] = (9*yPtoA*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=101
	expCoeffs[296] = (9*yPtoA*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=-11 t,u,v=201
	expCoeffs[297] = (9*yPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-11 t,u,v=301
	expCoeffs[298] = 0;
	// m,mp=-11 t,u,v=011
	expCoeffs[299] = (9*(1 + 2*expGamma*zPtoA*zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-11 t,u,v=111
	expCoeffs[300] = (9*(zPtoA + zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-11 t,u,v=211
	expCoeffs[301] = 9/(16.*std::pow(expGamma,4));
	// m,mp=-11 t,u,v=021
	expCoeffs[302] = 0;
	// m,mp=-11 t,u,v=121
	expCoeffs[303] = 0;
	// m,mp=-11 t,u,v=031
	expCoeffs[304] = 0;
	// m,mp=-11 t,u,v=002
	expCoeffs[305] = 0;
	// m,mp=-11 t,u,v=102
	expCoeffs[306] = 0;
	// m,mp=-11 t,u,v=202
	expCoeffs[307] = 0;
	// m,mp=-11 t,u,v=012
	expCoeffs[308] = 0;
	// m,mp=-11 t,u,v=112
	expCoeffs[309] = 0;
	// m,mp=-11 t,u,v=022
	expCoeffs[310] = 0;
	// m,mp=-11 t,u,v=003
	expCoeffs[311] = 0;
	// m,mp=-11 t,u,v=103
	expCoeffs[312] = 0;
	// m,mp=-11 t,u,v=013
	expCoeffs[313] = 0;
	// m,mp=-11 t,u,v=004
	expCoeffs[314] = 0;
	// m,mp=-12 t,u,v=000
	expCoeffs[315] = (-9*(yPtoB + expGamma*yPtoA*(-xPtoB + yPtoB)*(xPtoB + yPtoB))*zPtoA)/expGamma;
	// m,mp=-12 t,u,v=100
	expCoeffs[316] = (-9*(yPtoB + expGamma*yPtoA*(-xPtoB + yPtoB)*(xPtoB + yPtoB)))/(2.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=200
	expCoeffs[317] = 0;
	// m,mp=-12 t,u,v=300
	expCoeffs[318] = 0;
	// m,mp=-12 t,u,v=400
	expCoeffs[319] = 0;
	// m,mp=-12 t,u,v=010
	expCoeffs[320] = (9*(-1 + expGamma*std::pow(xPtoB,2) - expGamma*yPtoB*(2*yPtoA + yPtoB))*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=110
	expCoeffs[321] = (3*(-3 + 3*expGamma*std::pow(xPtoB,2) - 3*expGamma*yPtoB*(2*yPtoA + yPtoB)))/(4.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=210
	expCoeffs[322] = 0;
	// m,mp=-12 t,u,v=310
	expCoeffs[323] = 0;
	// m,mp=-12 t,u,v=020
	expCoeffs[324] = (-9*(yPtoA + 2*yPtoB)*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=120
	expCoeffs[325] = (-9*(yPtoA + 2*yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=220
	expCoeffs[326] = 0;
	// m,mp=-12 t,u,v=030
	expCoeffs[327] = (-9*zPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=130
	expCoeffs[328] = -9/(16.*std::pow(expGamma,4));
	// m,mp=-12 t,u,v=040
	expCoeffs[329] = 0;
	// m,mp=-12 t,u,v=001
	expCoeffs[330] = (9*xPtoB*yPtoA*zPtoA)/expGamma;
	// m,mp=-12 t,u,v=101
	expCoeffs[331] = (9*xPtoB*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=201
	expCoeffs[332] = 0;
	// m,mp=-12 t,u,v=301
	expCoeffs[333] = 0;
	// m,mp=-12 t,u,v=011
	expCoeffs[334] = (9*xPtoB*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=111
	expCoeffs[335] = (9*xPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=211
	expCoeffs[336] = 0;
	// m,mp=-12 t,u,v=021
	expCoeffs[337] = 0;
	// m,mp=-12 t,u,v=121
	expCoeffs[338] = 0;
	// m,mp=-12 t,u,v=031
	expCoeffs[339] = 0;
	// m,mp=-12 t,u,v=002
	expCoeffs[340] = (9*yPtoA*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=-12 t,u,v=102
	expCoeffs[341] = (9*yPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=202
	expCoeffs[342] = 0;
	// m,mp=-12 t,u,v=012
	expCoeffs[343] = (9*zPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=-12 t,u,v=112
	expCoeffs[344] = 9/(16.*std::pow(expGamma,4));
	// m,mp=-12 t,u,v=022
	expCoeffs[345] = 0;
	// m,mp=-12 t,u,v=003
	expCoeffs[346] = 0;
	// m,mp=-12 t,u,v=103
	expCoeffs[347] = 0;
	// m,mp=-12 t,u,v=013
	expCoeffs[348] = 0;
	// m,mp=-12 t,u,v=004
	expCoeffs[349] = 0;
	// m,mp=0-2 t,u,v=000
	expCoeffs[350] = (-3*(xPtoA*yPtoB + xPtoB*(yPtoA + expGamma*yPtoB*(d2PtoA - 3*std::pow(zPtoA,2)))))/expGamma;
	// m,mp=0-2 t,u,v=100
	expCoeffs[351] = (6*xPtoB*yPtoB*zPtoA)/expGamma;
	// m,mp=0-2 t,u,v=200
	expCoeffs[352] = (3*xPtoB*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=300
	expCoeffs[353] = 0;
	// m,mp=0-2 t,u,v=400
	expCoeffs[354] = 0;
	// m,mp=0-2 t,u,v=010
	expCoeffs[355] = (-3*(xPtoA + xPtoB + expGamma*xPtoB*(d2PtoA + 2*yPtoA*yPtoB - 3*std::pow(zPtoA,2))))/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=110
	expCoeffs[356] = (3*xPtoB*zPtoA)/std::pow(expGamma,2);
	// m,mp=0-2 t,u,v=210
	expCoeffs[357] = (3*xPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=310
	expCoeffs[358] = 0;
	// m,mp=0-2 t,u,v=020
	expCoeffs[359] = (-3*xPtoB*(2*yPtoA + yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=120
	expCoeffs[360] = 0;
	// m,mp=0-2 t,u,v=220
	expCoeffs[361] = 0;
	// m,mp=0-2 t,u,v=030
	expCoeffs[362] = (-3*xPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=130
	expCoeffs[363] = 0;
	// m,mp=0-2 t,u,v=040
	expCoeffs[364] = 0;
	// m,mp=0-2 t,u,v=001
	expCoeffs[365] = (-3*(yPtoA + yPtoB + expGamma*yPtoB*(d2PtoA + 2*xPtoA*xPtoB - 3*std::pow(zPtoA,2))))/(2.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=101
	expCoeffs[366] = (3*yPtoB*zPtoA)/std::pow(expGamma,2);
	// m,mp=0-2 t,u,v=201
	expCoeffs[367] = (3*yPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=301
	expCoeffs[368] = 0;
	// m,mp=0-2 t,u,v=011
	expCoeffs[369] = (-3*(2 + expGamma*(d2PtoA + 2*xPtoA*xPtoB + 2*yPtoA*yPtoB - 3*std::pow(zPtoA,2))))/(4.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=111
	expCoeffs[370] = (3*zPtoA)/(2.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=211
	expCoeffs[371] = 3/(8.*std::pow(expGamma,4));
	// m,mp=0-2 t,u,v=021
	expCoeffs[372] = (-3*(2*yPtoA + yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=121
	expCoeffs[373] = 0;
	// m,mp=0-2 t,u,v=031
	expCoeffs[374] = -3/(16.*std::pow(expGamma,4));
	// m,mp=0-2 t,u,v=002
	expCoeffs[375] = (-3*(2*xPtoA + xPtoB)*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=0-2 t,u,v=102
	expCoeffs[376] = 0;
	// m,mp=0-2 t,u,v=202
	expCoeffs[377] = 0;
	// m,mp=0-2 t,u,v=012
	expCoeffs[378] = (-3*(2*xPtoA + xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=112
	expCoeffs[379] = 0;
	// m,mp=0-2 t,u,v=022
	expCoeffs[380] = 0;
	// m,mp=0-2 t,u,v=003
	expCoeffs[381] = (-3*yPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=0-2 t,u,v=103
	expCoeffs[382] = 0;
	// m,mp=0-2 t,u,v=013
	expCoeffs[383] = -3/(16.*std::pow(expGamma,4));
	// m,mp=0-2 t,u,v=004
	expCoeffs[384] = 0;
	// m,mp=0-1 t,u,v=000
	expCoeffs[385] = (-3*(yPtoA*zPtoB + d2PtoA*expGamma*yPtoB*zPtoB - yPtoB*zPtoA*(2 + 3*expGamma*zPtoA*zPtoB)))/(2.*expGamma);
	// m,mp=0-1 t,u,v=100
	expCoeffs[386] = (-3*(yPtoA + yPtoB*(-2 + d2PtoA*expGamma - expGamma*zPtoA*(3*zPtoA + 4*zPtoB))))/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=200
	expCoeffs[387] = (3*yPtoB*(2*zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=300
	expCoeffs[388] = (3*yPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=400
	expCoeffs[389] = 0;
	// m,mp=0-1 t,u,v=010
	expCoeffs[390] = (-3*(-2*zPtoA + zPtoB + expGamma*(d2PtoA + 2*yPtoA*yPtoB - 3*std::pow(zPtoA,2))*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=110
	expCoeffs[391] = (3 - 3*expGamma*(d2PtoA + 2*yPtoA*yPtoB - zPtoA*(3*zPtoA + 4*zPtoB)))/(8.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=210
	expCoeffs[392] = (3*(2*zPtoA + zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=310
	expCoeffs[393] = 3/(16.*std::pow(expGamma,4));
	// m,mp=0-1 t,u,v=020
	expCoeffs[394] = (-3*(2*yPtoA + yPtoB)*zPtoB)/(8.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=120
	expCoeffs[395] = (-3*(2*yPtoA + yPtoB))/(16.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=220
	expCoeffs[396] = 0;
	// m,mp=0-1 t,u,v=030
	expCoeffs[397] = (-3*zPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=130
	expCoeffs[398] = -3/(32.*std::pow(expGamma,4));
	// m,mp=0-1 t,u,v=040
	expCoeffs[399] = 0;
	// m,mp=0-1 t,u,v=001
	expCoeffs[400] = (-3*xPtoA*yPtoB*zPtoB)/(2.*expGamma);
	// m,mp=0-1 t,u,v=101
	expCoeffs[401] = (-3*xPtoA*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=201
	expCoeffs[402] = 0;
	// m,mp=0-1 t,u,v=301
	expCoeffs[403] = 0;
	// m,mp=0-1 t,u,v=011
	expCoeffs[404] = (-3*xPtoA*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=111
	expCoeffs[405] = (-3*xPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=211
	expCoeffs[406] = 0;
	// m,mp=0-1 t,u,v=021
	expCoeffs[407] = 0;
	// m,mp=0-1 t,u,v=121
	expCoeffs[408] = 0;
	// m,mp=0-1 t,u,v=031
	expCoeffs[409] = 0;
	// m,mp=0-1 t,u,v=002
	expCoeffs[410] = (-3*yPtoB*zPtoB)/(8.*std::pow(expGamma,2));
	// m,mp=0-1 t,u,v=102
	expCoeffs[411] = (-3*yPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=202
	expCoeffs[412] = 0;
	// m,mp=0-1 t,u,v=012
	expCoeffs[413] = (-3*zPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=0-1 t,u,v=112
	expCoeffs[414] = -3/(32.*std::pow(expGamma,4));
	// m,mp=0-1 t,u,v=022
	expCoeffs[415] = 0;
	// m,mp=0-1 t,u,v=003
	expCoeffs[416] = 0;
	// m,mp=0-1 t,u,v=103
	expCoeffs[417] = 0;
	// m,mp=0-1 t,u,v=013
	expCoeffs[418] = 0;
	// m,mp=0-1 t,u,v=004
	expCoeffs[419] = 0;
	// m,mp=00 t,u,v=000
	expCoeffs[420] = (3 + expGamma*(2*xPtoA*xPtoB + 2*yPtoA*yPtoB + 8*zPtoA*zPtoB + expGamma*(d2PtoA - 3*std::pow(zPtoA,2))*(d2PtoB - 3*std::pow(zPtoB,2))))/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=100
	expCoeffs[421] = (2*(zPtoA + zPtoB) - expGamma*(d2PtoB*zPtoA + zPtoB*(d2PtoA - 3*zPtoA*(zPtoA + zPtoB))))/(2.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=200
	expCoeffs[422] = (4 - expGamma*(d2PtoA + d2PtoB - 3*std::pow(zPtoA,2) - 8*zPtoA*zPtoB - 3*std::pow(zPtoB,2)))/(8.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=300
	expCoeffs[423] = (zPtoA + zPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=400
	expCoeffs[424] = 1/(16.*std::pow(expGamma,4));
	// m,mp=00 t,u,v=010
	expCoeffs[425] = (yPtoA + yPtoB + expGamma*yPtoB*(d2PtoA - 3*std::pow(zPtoA,2)) + expGamma*yPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=110
	expCoeffs[426] = -0.5*(yPtoB*zPtoA + yPtoA*zPtoB)/std::pow(expGamma,2);
	// m,mp=00 t,u,v=210
	expCoeffs[427] = -0.125*(yPtoA + yPtoB)/std::pow(expGamma,3);
	// m,mp=00 t,u,v=310
	expCoeffs[428] = 0;
	// m,mp=00 t,u,v=020
	expCoeffs[429] = (2 + expGamma*(d2PtoA + d2PtoB + 4*yPtoA*yPtoB - 3*(std::pow(zPtoA,2) + std::pow(zPtoB,2))))/(16.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=120
	expCoeffs[430] = -0.125*(zPtoA + zPtoB)/std::pow(expGamma,3);
	// m,mp=00 t,u,v=220
	expCoeffs[431] = -0.0625*1/std::pow(expGamma,4);
	// m,mp=00 t,u,v=030
	expCoeffs[432] = (yPtoA + yPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=130
	expCoeffs[433] = 0;
	// m,mp=00 t,u,v=040
	expCoeffs[434] = 1/(64.*std::pow(expGamma,4));
	// m,mp=00 t,u,v=001
	expCoeffs[435] = (xPtoA + xPtoB + expGamma*xPtoB*(d2PtoA - 3*std::pow(zPtoA,2)) + expGamma*xPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=101
	expCoeffs[436] = -0.5*(xPtoB*zPtoA + xPtoA*zPtoB)/std::pow(expGamma,2);
	// m,mp=00 t,u,v=201
	expCoeffs[437] = -0.125*(xPtoA + xPtoB)/std::pow(expGamma,3);
	// m,mp=00 t,u,v=301
	expCoeffs[438] = 0;
	// m,mp=00 t,u,v=011
	expCoeffs[439] = (xPtoB*yPtoA + xPtoA*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=00 t,u,v=111
	expCoeffs[440] = 0;
	// m,mp=00 t,u,v=211
	expCoeffs[441] = 0;
	// m,mp=00 t,u,v=021
	expCoeffs[442] = (xPtoA + xPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=121
	expCoeffs[443] = 0;
	// m,mp=00 t,u,v=031
	expCoeffs[444] = 0;
	// m,mp=00 t,u,v=002
	expCoeffs[445] = (2 + expGamma*(d2PtoA + d2PtoB + 4*xPtoA*xPtoB - 3*(std::pow(zPtoA,2) + std::pow(zPtoB,2))))/(16.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=102
	expCoeffs[446] = -0.125*(zPtoA + zPtoB)/std::pow(expGamma,3);
	// m,mp=00 t,u,v=202
	expCoeffs[447] = -0.0625*1/std::pow(expGamma,4);
	// m,mp=00 t,u,v=012
	expCoeffs[448] = (yPtoA + yPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=112
	expCoeffs[449] = 0;
	// m,mp=00 t,u,v=022
	expCoeffs[450] = 1/(32.*std::pow(expGamma,4));
	// m,mp=00 t,u,v=003
	expCoeffs[451] = (xPtoA + xPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=00 t,u,v=103
	expCoeffs[452] = 0;
	// m,mp=00 t,u,v=013
	expCoeffs[453] = 0;
	// m,mp=00 t,u,v=004
	expCoeffs[454] = 1/(64.*std::pow(expGamma,4));
	// m,mp=01 t,u,v=000
	expCoeffs[455] = (-3*(xPtoA*zPtoB + d2PtoA*expGamma*xPtoB*zPtoB - xPtoB*zPtoA*(2 + 3*expGamma*zPtoA*zPtoB)))/(2.*expGamma);
	// m,mp=01 t,u,v=100
	expCoeffs[456] = (-3*(xPtoA + xPtoB*(-2 + d2PtoA*expGamma - expGamma*zPtoA*(3*zPtoA + 4*zPtoB))))/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=200
	expCoeffs[457] = (3*xPtoB*(2*zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=300
	expCoeffs[458] = (3*xPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=400
	expCoeffs[459] = 0;
	// m,mp=01 t,u,v=010
	expCoeffs[460] = (-3*xPtoB*yPtoA*zPtoB)/(2.*expGamma);
	// m,mp=01 t,u,v=110
	expCoeffs[461] = (-3*xPtoB*yPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=210
	expCoeffs[462] = 0;
	// m,mp=01 t,u,v=310
	expCoeffs[463] = 0;
	// m,mp=01 t,u,v=020
	expCoeffs[464] = (-3*xPtoB*zPtoB)/(8.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=120
	expCoeffs[465] = (-3*xPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=220
	expCoeffs[466] = 0;
	// m,mp=01 t,u,v=030
	expCoeffs[467] = 0;
	// m,mp=01 t,u,v=130
	expCoeffs[468] = 0;
	// m,mp=01 t,u,v=040
	expCoeffs[469] = 0;
	// m,mp=01 t,u,v=001
	expCoeffs[470] = (-3*(-2*zPtoA + zPtoB + expGamma*(d2PtoA + 2*xPtoA*xPtoB - 3*std::pow(zPtoA,2))*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=101
	expCoeffs[471] = (3 - 3*expGamma*(d2PtoA + 2*xPtoA*xPtoB - zPtoA*(3*zPtoA + 4*zPtoB)))/(8.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=201
	expCoeffs[472] = (3*(2*zPtoA + zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=301
	expCoeffs[473] = 3/(16.*std::pow(expGamma,4));
	// m,mp=01 t,u,v=011
	expCoeffs[474] = (-3*yPtoA*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=111
	expCoeffs[475] = (-3*yPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=211
	expCoeffs[476] = 0;
	// m,mp=01 t,u,v=021
	expCoeffs[477] = (-3*zPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=121
	expCoeffs[478] = -3/(32.*std::pow(expGamma,4));
	// m,mp=01 t,u,v=031
	expCoeffs[479] = 0;
	// m,mp=01 t,u,v=002
	expCoeffs[480] = (-3*(2*xPtoA + xPtoB)*zPtoB)/(8.*std::pow(expGamma,2));
	// m,mp=01 t,u,v=102
	expCoeffs[481] = (-3*(2*xPtoA + xPtoB))/(16.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=202
	expCoeffs[482] = 0;
	// m,mp=01 t,u,v=012
	expCoeffs[483] = 0;
	// m,mp=01 t,u,v=112
	expCoeffs[484] = 0;
	// m,mp=01 t,u,v=022
	expCoeffs[485] = 0;
	// m,mp=01 t,u,v=003
	expCoeffs[486] = (-3*zPtoB)/(16.*std::pow(expGamma,3));
	// m,mp=01 t,u,v=103
	expCoeffs[487] = -3/(32.*std::pow(expGamma,4));
	// m,mp=01 t,u,v=013
	expCoeffs[488] = 0;
	// m,mp=01 t,u,v=004
	expCoeffs[489] = 0;
	// m,mp=02 t,u,v=000
	expCoeffs[490] = (-6*xPtoA*xPtoB + 6*yPtoA*yPtoB + 3*expGamma*(-xPtoB + yPtoB)*(xPtoB + yPtoB)*(d2PtoA - 3*std::pow(zPtoA,2)))/(2.*expGamma);
	// m,mp=02 t,u,v=100
	expCoeffs[491] = (3*(xPtoB - yPtoB)*(xPtoB + yPtoB)*zPtoA)/expGamma;
	// m,mp=02 t,u,v=200
	expCoeffs[492] = (3*(xPtoB - yPtoB)*(xPtoB + yPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=300
	expCoeffs[493] = 0;
	// m,mp=02 t,u,v=400
	expCoeffs[494] = 0;
	// m,mp=02 t,u,v=010
	expCoeffs[495] = (3*(yPtoA - expGamma*std::pow(xPtoB,2)*yPtoA + yPtoB + expGamma*yPtoB*(d2PtoA + yPtoA*yPtoB - 3*std::pow(zPtoA,2))))/(2.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=110
	expCoeffs[496] = (-3*yPtoB*zPtoA)/std::pow(expGamma,2);
	// m,mp=02 t,u,v=210
	expCoeffs[497] = (-3*yPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=310
	expCoeffs[498] = 0;
	// m,mp=02 t,u,v=020
	expCoeffs[499] = (6 + 3*expGamma*(d2PtoA - std::pow(xPtoB,2) + 4*yPtoA*yPtoB + std::pow(yPtoB,2) - 3*std::pow(zPtoA,2)))/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=120
	expCoeffs[500] = (-3*zPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=220
	expCoeffs[501] = -3/(16.*std::pow(expGamma,4));
	// m,mp=02 t,u,v=030
	expCoeffs[502] = (3*(yPtoA + yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=130
	expCoeffs[503] = 0;
	// m,mp=02 t,u,v=040
	expCoeffs[504] = 3/(32.*std::pow(expGamma,4));
	// m,mp=02 t,u,v=001
	expCoeffs[505] = (-3*(xPtoA + xPtoB + expGamma*xPtoA*(xPtoB - yPtoB)*(xPtoB + yPtoB) + expGamma*xPtoB*(d2PtoA - 3*std::pow(zPtoA,2))))/(2.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=101
	expCoeffs[506] = (3*xPtoB*zPtoA)/std::pow(expGamma,2);
	// m,mp=02 t,u,v=201
	expCoeffs[507] = (3*xPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=301
	expCoeffs[508] = 0;
	// m,mp=02 t,u,v=011
	expCoeffs[509] = (3*(-(xPtoB*yPtoA) + xPtoA*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=02 t,u,v=111
	expCoeffs[510] = 0;
	// m,mp=02 t,u,v=211
	expCoeffs[511] = 0;
	// m,mp=02 t,u,v=021
	expCoeffs[512] = (3*(xPtoA - xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=121
	expCoeffs[513] = 0;
	// m,mp=02 t,u,v=031
	expCoeffs[514] = 0;
	// m,mp=02 t,u,v=002
	expCoeffs[515] = -0.125*(6 + 3*expGamma*(d2PtoA + 4*xPtoA*xPtoB + std::pow(xPtoB,2) - std::pow(yPtoB,2) - 3*std::pow(zPtoA,2)))/std::pow(expGamma,3);
	// m,mp=02 t,u,v=102
	expCoeffs[516] = (3*zPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=202
	expCoeffs[517] = 3/(16.*std::pow(expGamma,4));
	// m,mp=02 t,u,v=012
	expCoeffs[518] = (-3*(yPtoA - yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=112
	expCoeffs[519] = 0;
	// m,mp=02 t,u,v=022
	expCoeffs[520] = 0;
	// m,mp=02 t,u,v=003
	expCoeffs[521] = (-3*(xPtoA + xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=02 t,u,v=103
	expCoeffs[522] = 0;
	// m,mp=02 t,u,v=013
	expCoeffs[523] = 0;
	// m,mp=02 t,u,v=004
	expCoeffs[524] = -3/(32.*std::pow(expGamma,4));
	// m,mp=1-2 t,u,v=000
	expCoeffs[525] = 9*(1/expGamma + 2*xPtoA*xPtoB)*yPtoB*zPtoA;
	// m,mp=1-2 t,u,v=100
	expCoeffs[526] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=200
	expCoeffs[527] = 0;
	// m,mp=1-2 t,u,v=300
	expCoeffs[528] = 0;
	// m,mp=1-2 t,u,v=400
	expCoeffs[529] = 0;
	// m,mp=1-2 t,u,v=010
	expCoeffs[530] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=110
	expCoeffs[531] = (9*(1 + 2*expGamma*xPtoA*xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=1-2 t,u,v=210
	expCoeffs[532] = 0;
	// m,mp=1-2 t,u,v=310
	expCoeffs[533] = 0;
	// m,mp=1-2 t,u,v=020
	expCoeffs[534] = 0;
	// m,mp=1-2 t,u,v=120
	expCoeffs[535] = 0;
	// m,mp=1-2 t,u,v=220
	expCoeffs[536] = 0;
	// m,mp=1-2 t,u,v=030
	expCoeffs[537] = 0;
	// m,mp=1-2 t,u,v=130
	expCoeffs[538] = 0;
	// m,mp=1-2 t,u,v=040
	expCoeffs[539] = 0;
	// m,mp=1-2 t,u,v=001
	expCoeffs[540] = (9*(xPtoA + xPtoB)*yPtoB*zPtoA)/expGamma;
	// m,mp=1-2 t,u,v=101
	expCoeffs[541] = (9*(xPtoA + xPtoB)*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=201
	expCoeffs[542] = 0;
	// m,mp=1-2 t,u,v=301
	expCoeffs[543] = 0;
	// m,mp=1-2 t,u,v=011
	expCoeffs[544] = (9*(xPtoA + xPtoB)*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=111
	expCoeffs[545] = (9*(xPtoA + xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=1-2 t,u,v=211
	expCoeffs[546] = 0;
	// m,mp=1-2 t,u,v=021
	expCoeffs[547] = 0;
	// m,mp=1-2 t,u,v=121
	expCoeffs[548] = 0;
	// m,mp=1-2 t,u,v=031
	expCoeffs[549] = 0;
	// m,mp=1-2 t,u,v=002
	expCoeffs[550] = (9*yPtoB*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=1-2 t,u,v=102
	expCoeffs[551] = (9*yPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=1-2 t,u,v=202
	expCoeffs[552] = 0;
	// m,mp=1-2 t,u,v=012
	expCoeffs[553] = (9*zPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=1-2 t,u,v=112
	expCoeffs[554] = 9/(8.*std::pow(expGamma,4));
	// m,mp=1-2 t,u,v=022
	expCoeffs[555] = 0;
	// m,mp=1-2 t,u,v=003
	expCoeffs[556] = 0;
	// m,mp=1-2 t,u,v=103
	expCoeffs[557] = 0;
	// m,mp=1-2 t,u,v=013
	expCoeffs[558] = 0;
	// m,mp=1-2 t,u,v=004
	expCoeffs[559] = 0;
	// m,mp=1-1 t,u,v=000
	expCoeffs[560] = (9*xPtoA*yPtoB*(1/expGamma + 2*zPtoA*zPtoB))/2.;
	// m,mp=1-1 t,u,v=100
	expCoeffs[561] = (9*xPtoA*yPtoB*(zPtoA + zPtoB))/(2.*expGamma);
	// m,mp=1-1 t,u,v=200
	expCoeffs[562] = (9*xPtoA*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=300
	expCoeffs[563] = 0;
	// m,mp=1-1 t,u,v=400
	expCoeffs[564] = 0;
	// m,mp=1-1 t,u,v=010
	expCoeffs[565] = (9*xPtoA*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=110
	expCoeffs[566] = (9*xPtoA*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=210
	expCoeffs[567] = (9*xPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=1-1 t,u,v=310
	expCoeffs[568] = 0;
	// m,mp=1-1 t,u,v=020
	expCoeffs[569] = 0;
	// m,mp=1-1 t,u,v=120
	expCoeffs[570] = 0;
	// m,mp=1-1 t,u,v=220
	expCoeffs[571] = 0;
	// m,mp=1-1 t,u,v=030
	expCoeffs[572] = 0;
	// m,mp=1-1 t,u,v=130
	expCoeffs[573] = 0;
	// m,mp=1-1 t,u,v=040
	expCoeffs[574] = 0;
	// m,mp=1-1 t,u,v=001
	expCoeffs[575] = (9*yPtoB*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=101
	expCoeffs[576] = (9*yPtoB*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=1-1 t,u,v=201
	expCoeffs[577] = (9*yPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=1-1 t,u,v=301
	expCoeffs[578] = 0;
	// m,mp=1-1 t,u,v=011
	expCoeffs[579] = (9*(1 + 2*expGamma*zPtoA*zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=1-1 t,u,v=111
	expCoeffs[580] = (9*(zPtoA + zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=1-1 t,u,v=211
	expCoeffs[581] = 9/(16.*std::pow(expGamma,4));
	// m,mp=1-1 t,u,v=021
	expCoeffs[582] = 0;
	// m,mp=1-1 t,u,v=121
	expCoeffs[583] = 0;
	// m,mp=1-1 t,u,v=031
	expCoeffs[584] = 0;
	// m,mp=1-1 t,u,v=002
	expCoeffs[585] = 0;
	// m,mp=1-1 t,u,v=102
	expCoeffs[586] = 0;
	// m,mp=1-1 t,u,v=202
	expCoeffs[587] = 0;
	// m,mp=1-1 t,u,v=012
	expCoeffs[588] = 0;
	// m,mp=1-1 t,u,v=112
	expCoeffs[589] = 0;
	// m,mp=1-1 t,u,v=022
	expCoeffs[590] = 0;
	// m,mp=1-1 t,u,v=003
	expCoeffs[591] = 0;
	// m,mp=1-1 t,u,v=103
	expCoeffs[592] = 0;
	// m,mp=1-1 t,u,v=013
	expCoeffs[593] = 0;
	// m,mp=1-1 t,u,v=004
	expCoeffs[594] = 0;
	// m,mp=10 t,u,v=000
	expCoeffs[595] = (3*(2*xPtoA*zPtoB - zPtoA*(xPtoB + expGamma*xPtoA*(d2PtoB - 3*std::pow(zPtoB,2)))))/(2.*expGamma);
	// m,mp=10 t,u,v=100
	expCoeffs[596] = (-3*(xPtoB + xPtoA*(-2 + d2PtoB*expGamma - expGamma*zPtoB*(4*zPtoA + 3*zPtoB))))/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=200
	expCoeffs[597] = (3*xPtoA*(zPtoA + 2*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=300
	expCoeffs[598] = (3*xPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=400
	expCoeffs[599] = 0;
	// m,mp=10 t,u,v=010
	expCoeffs[600] = (-3*xPtoA*yPtoB*zPtoA)/(2.*expGamma);
	// m,mp=10 t,u,v=110
	expCoeffs[601] = (-3*xPtoA*yPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=210
	expCoeffs[602] = 0;
	// m,mp=10 t,u,v=310
	expCoeffs[603] = 0;
	// m,mp=10 t,u,v=020
	expCoeffs[604] = (-3*xPtoA*zPtoA)/(8.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=120
	expCoeffs[605] = (-3*xPtoA)/(16.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=220
	expCoeffs[606] = 0;
	// m,mp=10 t,u,v=030
	expCoeffs[607] = 0;
	// m,mp=10 t,u,v=130
	expCoeffs[608] = 0;
	// m,mp=10 t,u,v=040
	expCoeffs[609] = 0;
	// m,mp=10 t,u,v=001
	expCoeffs[610] = (-3*(zPtoA - 2*zPtoB + expGamma*zPtoA*(d2PtoB + 2*xPtoA*xPtoB - 3*std::pow(zPtoB,2))))/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=101
	expCoeffs[611] = (3 - 3*expGamma*(d2PtoB + 2*xPtoA*xPtoB - zPtoB*(4*zPtoA + 3*zPtoB)))/(8.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=201
	expCoeffs[612] = (3*(zPtoA + 2*zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=301
	expCoeffs[613] = 3/(16.*std::pow(expGamma,4));
	// m,mp=10 t,u,v=011
	expCoeffs[614] = (-3*yPtoB*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=111
	expCoeffs[615] = (-3*yPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=211
	expCoeffs[616] = 0;
	// m,mp=10 t,u,v=021
	expCoeffs[617] = (-3*zPtoA)/(16.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=121
	expCoeffs[618] = -3/(32.*std::pow(expGamma,4));
	// m,mp=10 t,u,v=031
	expCoeffs[619] = 0;
	// m,mp=10 t,u,v=002
	expCoeffs[620] = (-3*(xPtoA + 2*xPtoB)*zPtoA)/(8.*std::pow(expGamma,2));
	// m,mp=10 t,u,v=102
	expCoeffs[621] = (-3*(xPtoA + 2*xPtoB))/(16.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=202
	expCoeffs[622] = 0;
	// m,mp=10 t,u,v=012
	expCoeffs[623] = 0;
	// m,mp=10 t,u,v=112
	expCoeffs[624] = 0;
	// m,mp=10 t,u,v=022
	expCoeffs[625] = 0;
	// m,mp=10 t,u,v=003
	expCoeffs[626] = (-3*zPtoA)/(16.*std::pow(expGamma,3));
	// m,mp=10 t,u,v=103
	expCoeffs[627] = -3/(32.*std::pow(expGamma,4));
	// m,mp=10 t,u,v=013
	expCoeffs[628] = 0;
	// m,mp=10 t,u,v=004
	expCoeffs[629] = 0;
	// m,mp=11 t,u,v=000
	expCoeffs[630] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=100
	expCoeffs[631] = (9*(1 + 2*expGamma*xPtoA*xPtoB)*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=200
	expCoeffs[632] = (9*(1 + 2*expGamma*xPtoA*xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=11 t,u,v=300
	expCoeffs[633] = 0;
	// m,mp=11 t,u,v=400
	expCoeffs[634] = 0;
	// m,mp=11 t,u,v=010
	expCoeffs[635] = 0;
	// m,mp=11 t,u,v=110
	expCoeffs[636] = 0;
	// m,mp=11 t,u,v=210
	expCoeffs[637] = 0;
	// m,mp=11 t,u,v=310
	expCoeffs[638] = 0;
	// m,mp=11 t,u,v=020
	expCoeffs[639] = 0;
	// m,mp=11 t,u,v=120
	expCoeffs[640] = 0;
	// m,mp=11 t,u,v=220
	expCoeffs[641] = 0;
	// m,mp=11 t,u,v=030
	expCoeffs[642] = 0;
	// m,mp=11 t,u,v=130
	expCoeffs[643] = 0;
	// m,mp=11 t,u,v=040
	expCoeffs[644] = 0;
	// m,mp=11 t,u,v=001
	expCoeffs[645] = (9*(xPtoA + xPtoB)*(1 + 2*expGamma*zPtoA*zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=101
	expCoeffs[646] = (9*(xPtoA + xPtoB)*(zPtoA + zPtoB))/(4.*std::pow(expGamma,2));
	// m,mp=11 t,u,v=201
	expCoeffs[647] = (9*(xPtoA + xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=11 t,u,v=301
	expCoeffs[648] = 0;
	// m,mp=11 t,u,v=011
	expCoeffs[649] = 0;
	// m,mp=11 t,u,v=111
	expCoeffs[650] = 0;
	// m,mp=11 t,u,v=211
	expCoeffs[651] = 0;
	// m,mp=11 t,u,v=021
	expCoeffs[652] = 0;
	// m,mp=11 t,u,v=121
	expCoeffs[653] = 0;
	// m,mp=11 t,u,v=031
	expCoeffs[654] = 0;
	// m,mp=11 t,u,v=002
	expCoeffs[655] = (9*(1 + 2*expGamma*zPtoA*zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=11 t,u,v=102
	expCoeffs[656] = (9*(zPtoA + zPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=11 t,u,v=202
	expCoeffs[657] = 9/(16.*std::pow(expGamma,4));
	// m,mp=11 t,u,v=012
	expCoeffs[658] = 0;
	// m,mp=11 t,u,v=112
	expCoeffs[659] = 0;
	// m,mp=11 t,u,v=022
	expCoeffs[660] = 0;
	// m,mp=11 t,u,v=003
	expCoeffs[661] = 0;
	// m,mp=11 t,u,v=103
	expCoeffs[662] = 0;
	// m,mp=11 t,u,v=013
	expCoeffs[663] = 0;
	// m,mp=11 t,u,v=004
	expCoeffs[664] = 0;
	// m,mp=12 t,u,v=000
	expCoeffs[665] = (9*(xPtoB + expGamma*xPtoA*(xPtoB - yPtoB)*(xPtoB + yPtoB))*zPtoA)/expGamma;
	// m,mp=12 t,u,v=100
	expCoeffs[666] = (9*(xPtoB + expGamma*xPtoA*(xPtoB - yPtoB)*(xPtoB + yPtoB)))/(2.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=200
	expCoeffs[667] = 0;
	// m,mp=12 t,u,v=300
	expCoeffs[668] = 0;
	// m,mp=12 t,u,v=400
	expCoeffs[669] = 0;
	// m,mp=12 t,u,v=010
	expCoeffs[670] = (-9*xPtoA*yPtoB*zPtoA)/expGamma;
	// m,mp=12 t,u,v=110
	expCoeffs[671] = (-9*xPtoA*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=210
	expCoeffs[672] = 0;
	// m,mp=12 t,u,v=310
	expCoeffs[673] = 0;
	// m,mp=12 t,u,v=020
	expCoeffs[674] = (-9*xPtoA*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=120
	expCoeffs[675] = (-9*xPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=220
	expCoeffs[676] = 0;
	// m,mp=12 t,u,v=030
	expCoeffs[677] = 0;
	// m,mp=12 t,u,v=130
	expCoeffs[678] = 0;
	// m,mp=12 t,u,v=040
	expCoeffs[679] = 0;
	// m,mp=12 t,u,v=001
	expCoeffs[680] = (9*(1 + 2*expGamma*xPtoA*xPtoB + expGamma*(xPtoB - yPtoB)*(xPtoB + yPtoB))*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=101
	expCoeffs[681] = (9*(1 + 2*expGamma*xPtoA*xPtoB + expGamma*(xPtoB - yPtoB)*(xPtoB + yPtoB)))/(4.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=201
	expCoeffs[682] = 0;
	// m,mp=12 t,u,v=301
	expCoeffs[683] = 0;
	// m,mp=12 t,u,v=011
	expCoeffs[684] = (-9*yPtoB*zPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=111
	expCoeffs[685] = (-9*yPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=211
	expCoeffs[686] = 0;
	// m,mp=12 t,u,v=021
	expCoeffs[687] = (-9*zPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=121
	expCoeffs[688] = -9/(16.*std::pow(expGamma,4));
	// m,mp=12 t,u,v=031
	expCoeffs[689] = 0;
	// m,mp=12 t,u,v=002
	expCoeffs[690] = (9*(xPtoA + 2*xPtoB)*zPtoA)/(4.*std::pow(expGamma,2));
	// m,mp=12 t,u,v=102
	expCoeffs[691] = (9*(xPtoA + 2*xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=202
	expCoeffs[692] = 0;
	// m,mp=12 t,u,v=012
	expCoeffs[693] = 0;
	// m,mp=12 t,u,v=112
	expCoeffs[694] = 0;
	// m,mp=12 t,u,v=022
	expCoeffs[695] = 0;
	// m,mp=12 t,u,v=003
	expCoeffs[696] = (9*zPtoA)/(8.*std::pow(expGamma,3));
	// m,mp=12 t,u,v=103
	expCoeffs[697] = 9/(16.*std::pow(expGamma,4));
	// m,mp=12 t,u,v=013
	expCoeffs[698] = 0;
	// m,mp=12 t,u,v=004
	expCoeffs[699] = 0;
	// m,mp=2-2 t,u,v=000
	expCoeffs[700] = (18*(xPtoA*yPtoB - xPtoB*(yPtoA + expGamma*(-xPtoA + yPtoA)*(xPtoA + yPtoA)*yPtoB)))/expGamma;
	// m,mp=2-2 t,u,v=100
	expCoeffs[701] = 0;
	// m,mp=2-2 t,u,v=200
	expCoeffs[702] = 0;
	// m,mp=2-2 t,u,v=300
	expCoeffs[703] = 0;
	// m,mp=2-2 t,u,v=400
	expCoeffs[704] = 0;
	// m,mp=2-2 t,u,v=010
	expCoeffs[705] = (9*(xPtoA + expGamma*std::pow(xPtoA,2)*xPtoB - xPtoB*(1 + expGamma*yPtoA*(yPtoA + 2*yPtoB))))/std::pow(expGamma,2);
	// m,mp=2-2 t,u,v=110
	expCoeffs[706] = 0;
	// m,mp=2-2 t,u,v=210
	expCoeffs[707] = 0;
	// m,mp=2-2 t,u,v=310
	expCoeffs[708] = 0;
	// m,mp=2-2 t,u,v=020
	expCoeffs[709] = (-9*xPtoB*(2*yPtoA + yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=2-2 t,u,v=120
	expCoeffs[710] = 0;
	// m,mp=2-2 t,u,v=220
	expCoeffs[711] = 0;
	// m,mp=2-2 t,u,v=030
	expCoeffs[712] = (-9*xPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=2-2 t,u,v=130
	expCoeffs[713] = 0;
	// m,mp=2-2 t,u,v=040
	expCoeffs[714] = 0;
	// m,mp=2-2 t,u,v=001
	expCoeffs[715] = (9*(yPtoB + expGamma*xPtoA*(xPtoA + 2*xPtoB)*yPtoB - yPtoA*(1 + expGamma*yPtoA*yPtoB)))/std::pow(expGamma,2);
	// m,mp=2-2 t,u,v=101
	expCoeffs[716] = 0;
	// m,mp=2-2 t,u,v=201
	expCoeffs[717] = 0;
	// m,mp=2-2 t,u,v=301
	expCoeffs[718] = 0;
	// m,mp=2-2 t,u,v=011
	expCoeffs[719] = (9*(std::pow(xPtoA,2) + 2*xPtoA*xPtoB - yPtoA*(yPtoA + 2*yPtoB)))/(2.*std::pow(expGamma,2));
	// m,mp=2-2 t,u,v=111
	expCoeffs[720] = 0;
	// m,mp=2-2 t,u,v=211
	expCoeffs[721] = 0;
	// m,mp=2-2 t,u,v=021
	expCoeffs[722] = (-9*(2*yPtoA + yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=2-2 t,u,v=121
	expCoeffs[723] = 0;
	// m,mp=2-2 t,u,v=031
	expCoeffs[724] = -9/(8.*std::pow(expGamma,4));
	// m,mp=2-2 t,u,v=002
	expCoeffs[725] = (9*(2*xPtoA + xPtoB)*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=2-2 t,u,v=102
	expCoeffs[726] = 0;
	// m,mp=2-2 t,u,v=202
	expCoeffs[727] = 0;
	// m,mp=2-2 t,u,v=012
	expCoeffs[728] = (9*(2*xPtoA + xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=2-2 t,u,v=112
	expCoeffs[729] = 0;
	// m,mp=2-2 t,u,v=022
	expCoeffs[730] = 0;
	// m,mp=2-2 t,u,v=003
	expCoeffs[731] = (9*yPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=2-2 t,u,v=103
	expCoeffs[732] = 0;
	// m,mp=2-2 t,u,v=013
	expCoeffs[733] = 9/(8.*std::pow(expGamma,4));
	// m,mp=2-2 t,u,v=004
	expCoeffs[734] = 0;
	// m,mp=2-1 t,u,v=000
	expCoeffs[735] = (-9*(yPtoA + expGamma*(-xPtoA + yPtoA)*(xPtoA + yPtoA)*yPtoB)*zPtoB)/expGamma;
	// m,mp=2-1 t,u,v=100
	expCoeffs[736] = (-9*(yPtoA + expGamma*(-xPtoA + yPtoA)*(xPtoA + yPtoA)*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=200
	expCoeffs[737] = 0;
	// m,mp=2-1 t,u,v=300
	expCoeffs[738] = 0;
	// m,mp=2-1 t,u,v=400
	expCoeffs[739] = 0;
	// m,mp=2-1 t,u,v=010
	expCoeffs[740] = (9*(-1 + expGamma*std::pow(xPtoA,2) - expGamma*yPtoA*(yPtoA + 2*yPtoB))*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=110
	expCoeffs[741] = (-9 + 9*expGamma*std::pow(xPtoA,2) - 9*expGamma*yPtoA*(yPtoA + 2*yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=210
	expCoeffs[742] = 0;
	// m,mp=2-1 t,u,v=310
	expCoeffs[743] = 0;
	// m,mp=2-1 t,u,v=020
	expCoeffs[744] = (-9*(2*yPtoA + yPtoB)*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=120
	expCoeffs[745] = (-9*(2*yPtoA + yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=220
	expCoeffs[746] = 0;
	// m,mp=2-1 t,u,v=030
	expCoeffs[747] = (-9*zPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=130
	expCoeffs[748] = -9/(16.*std::pow(expGamma,4));
	// m,mp=2-1 t,u,v=040
	expCoeffs[749] = 0;
	// m,mp=2-1 t,u,v=001
	expCoeffs[750] = (9*xPtoA*yPtoB*zPtoB)/expGamma;
	// m,mp=2-1 t,u,v=101
	expCoeffs[751] = (9*xPtoA*yPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=201
	expCoeffs[752] = 0;
	// m,mp=2-1 t,u,v=301
	expCoeffs[753] = 0;
	// m,mp=2-1 t,u,v=011
	expCoeffs[754] = (9*xPtoA*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=111
	expCoeffs[755] = (9*xPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=211
	expCoeffs[756] = 0;
	// m,mp=2-1 t,u,v=021
	expCoeffs[757] = 0;
	// m,mp=2-1 t,u,v=121
	expCoeffs[758] = 0;
	// m,mp=2-1 t,u,v=031
	expCoeffs[759] = 0;
	// m,mp=2-1 t,u,v=002
	expCoeffs[760] = (9*yPtoB*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=2-1 t,u,v=102
	expCoeffs[761] = (9*yPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=202
	expCoeffs[762] = 0;
	// m,mp=2-1 t,u,v=012
	expCoeffs[763] = (9*zPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=2-1 t,u,v=112
	expCoeffs[764] = 9/(16.*std::pow(expGamma,4));
	// m,mp=2-1 t,u,v=022
	expCoeffs[765] = 0;
	// m,mp=2-1 t,u,v=003
	expCoeffs[766] = 0;
	// m,mp=2-1 t,u,v=103
	expCoeffs[767] = 0;
	// m,mp=2-1 t,u,v=013
	expCoeffs[768] = 0;
	// m,mp=2-1 t,u,v=004
	expCoeffs[769] = 0;
	// m,mp=20 t,u,v=000
	expCoeffs[770] = (-6*xPtoA*xPtoB + 6*yPtoA*yPtoB - 3*expGamma*(xPtoA - yPtoA)*(xPtoA + yPtoA)*(d2PtoB - 3*std::pow(zPtoB,2)))/(2.*expGamma);
	// m,mp=20 t,u,v=100
	expCoeffs[771] = (3*(xPtoA - yPtoA)*(xPtoA + yPtoA)*zPtoB)/expGamma;
	// m,mp=20 t,u,v=200
	expCoeffs[772] = (3*(xPtoA - yPtoA)*(xPtoA + yPtoA))/(4.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=300
	expCoeffs[773] = 0;
	// m,mp=20 t,u,v=400
	expCoeffs[774] = 0;
	// m,mp=20 t,u,v=010
	expCoeffs[775] = (3*(yPtoA + yPtoB - expGamma*std::pow(xPtoA,2)*yPtoB + expGamma*yPtoA*(d2PtoB + yPtoA*yPtoB - 3*std::pow(zPtoB,2))))/(2.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=110
	expCoeffs[776] = (-3*yPtoA*zPtoB)/std::pow(expGamma,2);
	// m,mp=20 t,u,v=210
	expCoeffs[777] = (-3*yPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=310
	expCoeffs[778] = 0;
	// m,mp=20 t,u,v=020
	expCoeffs[779] = (6 + 3*expGamma*(d2PtoB - std::pow(xPtoA,2) + std::pow(yPtoA,2) + 4*yPtoA*yPtoB - 3*std::pow(zPtoB,2)))/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=120
	expCoeffs[780] = (-3*zPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=220
	expCoeffs[781] = -3/(16.*std::pow(expGamma,4));
	// m,mp=20 t,u,v=030
	expCoeffs[782] = (3*(yPtoA + yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=130
	expCoeffs[783] = 0;
	// m,mp=20 t,u,v=040
	expCoeffs[784] = 3/(32.*std::pow(expGamma,4));
	// m,mp=20 t,u,v=001
	expCoeffs[785] = (-3*(xPtoA + xPtoB + expGamma*std::pow(xPtoA,2)*xPtoB - expGamma*xPtoB*std::pow(yPtoA,2) + expGamma*xPtoA*(d2PtoB - 3*std::pow(zPtoB,2))))/(2.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=101
	expCoeffs[786] = (3*xPtoA*zPtoB)/std::pow(expGamma,2);
	// m,mp=20 t,u,v=201
	expCoeffs[787] = (3*xPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=301
	expCoeffs[788] = 0;
	// m,mp=20 t,u,v=011
	expCoeffs[789] = (3*(xPtoB*yPtoA - xPtoA*yPtoB))/(2.*std::pow(expGamma,2));
	// m,mp=20 t,u,v=111
	expCoeffs[790] = 0;
	// m,mp=20 t,u,v=211
	expCoeffs[791] = 0;
	// m,mp=20 t,u,v=021
	expCoeffs[792] = (-3*(xPtoA - xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=121
	expCoeffs[793] = 0;
	// m,mp=20 t,u,v=031
	expCoeffs[794] = 0;
	// m,mp=20 t,u,v=002
	expCoeffs[795] = -0.125*(6 + 3*expGamma*(d2PtoB + std::pow(xPtoA,2) + 4*xPtoA*xPtoB - std::pow(yPtoA,2) - 3*std::pow(zPtoB,2)))/std::pow(expGamma,3);
	// m,mp=20 t,u,v=102
	expCoeffs[796] = (3*zPtoB)/(4.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=202
	expCoeffs[797] = 3/(16.*std::pow(expGamma,4));
	// m,mp=20 t,u,v=012
	expCoeffs[798] = (3*(yPtoA - yPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=112
	expCoeffs[799] = 0;
	// m,mp=20 t,u,v=022
	expCoeffs[800] = 0;
	// m,mp=20 t,u,v=003
	expCoeffs[801] = (-3*(xPtoA + xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=20 t,u,v=103
	expCoeffs[802] = 0;
	// m,mp=20 t,u,v=013
	expCoeffs[803] = 0;
	// m,mp=20 t,u,v=004
	expCoeffs[804] = -3/(32.*std::pow(expGamma,4));
	// m,mp=21 t,u,v=000
	expCoeffs[805] = (9*(xPtoA + expGamma*xPtoB*(xPtoA - yPtoA)*(xPtoA + yPtoA))*zPtoB)/expGamma;
	// m,mp=21 t,u,v=100
	expCoeffs[806] = (9*(xPtoA + expGamma*xPtoB*(xPtoA - yPtoA)*(xPtoA + yPtoA)))/(2.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=200
	expCoeffs[807] = 0;
	// m,mp=21 t,u,v=300
	expCoeffs[808] = 0;
	// m,mp=21 t,u,v=400
	expCoeffs[809] = 0;
	// m,mp=21 t,u,v=010
	expCoeffs[810] = (-9*xPtoB*yPtoA*zPtoB)/expGamma;
	// m,mp=21 t,u,v=110
	expCoeffs[811] = (-9*xPtoB*yPtoA)/(2.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=210
	expCoeffs[812] = 0;
	// m,mp=21 t,u,v=310
	expCoeffs[813] = 0;
	// m,mp=21 t,u,v=020
	expCoeffs[814] = (-9*xPtoB*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=120
	expCoeffs[815] = (-9*xPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=220
	expCoeffs[816] = 0;
	// m,mp=21 t,u,v=030
	expCoeffs[817] = 0;
	// m,mp=21 t,u,v=130
	expCoeffs[818] = 0;
	// m,mp=21 t,u,v=040
	expCoeffs[819] = 0;
	// m,mp=21 t,u,v=001
	expCoeffs[820] = (9*(1 + expGamma*(std::pow(xPtoA,2) + 2*xPtoA*xPtoB - std::pow(yPtoA,2)))*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=101
	expCoeffs[821] = (9*(1 + expGamma*(std::pow(xPtoA,2) + 2*xPtoA*xPtoB - std::pow(yPtoA,2))))/(4.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=201
	expCoeffs[822] = 0;
	// m,mp=21 t,u,v=301
	expCoeffs[823] = 0;
	// m,mp=21 t,u,v=011
	expCoeffs[824] = (-9*yPtoA*zPtoB)/(2.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=111
	expCoeffs[825] = (-9*yPtoA)/(4.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=211
	expCoeffs[826] = 0;
	// m,mp=21 t,u,v=021
	expCoeffs[827] = (-9*zPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=121
	expCoeffs[828] = -9/(16.*std::pow(expGamma,4));
	// m,mp=21 t,u,v=031
	expCoeffs[829] = 0;
	// m,mp=21 t,u,v=002
	expCoeffs[830] = (9*(2*xPtoA + xPtoB)*zPtoB)/(4.*std::pow(expGamma,2));
	// m,mp=21 t,u,v=102
	expCoeffs[831] = (9*(2*xPtoA + xPtoB))/(8.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=202
	expCoeffs[832] = 0;
	// m,mp=21 t,u,v=012
	expCoeffs[833] = 0;
	// m,mp=21 t,u,v=112
	expCoeffs[834] = 0;
	// m,mp=21 t,u,v=022
	expCoeffs[835] = 0;
	// m,mp=21 t,u,v=003
	expCoeffs[836] = (9*zPtoB)/(8.*std::pow(expGamma,3));
	// m,mp=21 t,u,v=103
	expCoeffs[837] = 9/(16.*std::pow(expGamma,4));
	// m,mp=21 t,u,v=013
	expCoeffs[838] = 0;
	// m,mp=21 t,u,v=004
	expCoeffs[839] = 0;
	// m,mp=22 t,u,v=000
	expCoeffs[840] = (9 + 9*expGamma*(2*xPtoA*xPtoB + 2*yPtoA*yPtoB + expGamma*(xPtoA - yPtoA)*(xPtoA + yPtoA)*(xPtoB - yPtoB)*(xPtoB + yPtoB)))/std::pow(expGamma,2);
	// m,mp=22 t,u,v=100
	expCoeffs[841] = 0;
	// m,mp=22 t,u,v=200
	expCoeffs[842] = 0;
	// m,mp=22 t,u,v=300
	expCoeffs[843] = 0;
	// m,mp=22 t,u,v=400
	expCoeffs[844] = 0;
	// m,mp=22 t,u,v=010
	expCoeffs[845] = (9*(yPtoA + yPtoB + expGamma*(-(std::pow(xPtoB,2)*yPtoA) - std::pow(xPtoA,2)*yPtoB + yPtoA*yPtoB*(yPtoA + yPtoB))))/std::pow(expGamma,2);
	// m,mp=22 t,u,v=110
	expCoeffs[846] = 0;
	// m,mp=22 t,u,v=210
	expCoeffs[847] = 0;
	// m,mp=22 t,u,v=310
	expCoeffs[848] = 0;
	// m,mp=22 t,u,v=020
	expCoeffs[849] = (18 + 9*expGamma*(-std::pow(xPtoA,2) - std::pow(xPtoB,2) + std::pow(yPtoA,2) + 4*yPtoA*yPtoB + std::pow(yPtoB,2)))/(4.*std::pow(expGamma,3));
	// m,mp=22 t,u,v=120
	expCoeffs[850] = 0;
	// m,mp=22 t,u,v=220
	expCoeffs[851] = 0;
	// m,mp=22 t,u,v=030
	expCoeffs[852] = (9*(yPtoA + yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=22 t,u,v=130
	expCoeffs[853] = 0;
	// m,mp=22 t,u,v=040
	expCoeffs[854] = 9/(16.*std::pow(expGamma,4));
	// m,mp=22 t,u,v=001
	expCoeffs[855] = (9*(xPtoA + xPtoB + expGamma*(xPtoA*xPtoB*(xPtoA + xPtoB) - xPtoB*std::pow(yPtoA,2) - xPtoA*std::pow(yPtoB,2))))/std::pow(expGamma,2);
	// m,mp=22 t,u,v=101
	expCoeffs[856] = 0;
	// m,mp=22 t,u,v=201
	expCoeffs[857] = 0;
	// m,mp=22 t,u,v=301
	expCoeffs[858] = 0;
	// m,mp=22 t,u,v=011
	expCoeffs[859] = (-9*(xPtoB*yPtoA + xPtoA*yPtoB))/std::pow(expGamma,2);
	// m,mp=22 t,u,v=111
	expCoeffs[860] = 0;
	// m,mp=22 t,u,v=211
	expCoeffs[861] = 0;
	// m,mp=22 t,u,v=021
	expCoeffs[862] = (-9*(xPtoA + xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=22 t,u,v=121
	expCoeffs[863] = 0;
	// m,mp=22 t,u,v=031
	expCoeffs[864] = 0;
	// m,mp=22 t,u,v=002
	expCoeffs[865] = (18 + 9*expGamma*(std::pow(xPtoA,2) + 4*xPtoA*xPtoB + std::pow(xPtoB,2) - std::pow(yPtoA,2) - std::pow(yPtoB,2)))/(4.*std::pow(expGamma,3));
	// m,mp=22 t,u,v=102
	expCoeffs[866] = 0;
	// m,mp=22 t,u,v=202
	expCoeffs[867] = 0;
	// m,mp=22 t,u,v=012
	expCoeffs[868] = (-9*(yPtoA + yPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=22 t,u,v=112
	expCoeffs[869] = 0;
	// m,mp=22 t,u,v=022
	expCoeffs[870] = -9/(8.*std::pow(expGamma,4));
	// m,mp=22 t,u,v=003
	expCoeffs[871] = (9*(xPtoA + xPtoB))/(4.*std::pow(expGamma,3));
	// m,mp=22 t,u,v=103
	expCoeffs[872] = 0;
	// m,mp=22 t,u,v=013
	expCoeffs[873] = 0;
	// m,mp=22 t,u,v=004
	expCoeffs[874] = 9/(16.*std::pow(expGamma,4));
	return expCoeffs;
}
