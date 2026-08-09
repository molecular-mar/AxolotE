#ifndef ECOEFMA_H
#define ECOEFMA_H

#include <array>

std::vector<double> ECoeff0And0(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff0And1(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff0And2(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff1And0(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff1And1(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff1And2(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff2And0(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff2And1(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
std::vector<double> ECoeff2And2(double expAlpha, double expBeta,         
	std::array<double,3> muCoords, std::array<double,3> nuCoords);
#endif