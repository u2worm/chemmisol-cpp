#include "chemical.h"
#include "newton.h"

#include <iostream>
#include <limits>

namespace newton {
	template<>
		double abs(const X& x) {
			double a = std::numeric_limits<double>::min();
			for(auto& v : x) {
				double abs_v = std::abs(v);
				if (abs_v > a)
					a = abs_v;
			}
			return a;
		}
}

const double ChemicalSystem::bulk_density = 1.17*gram/cm3;
const double ChemicalSystem::V = std::pow(1*cm, 3);
const double ChemicalSystem::mineral_weight  = V*bulk_density;

double deltaH(const X& x) {
	return x[0];
}
double dH(const X& x) {
	return x[1];
}
double dP(const X& x) {
	return x[2];
}
double dC(const X& x) {
	return x[3];
}

float trace(const M& m) {
	float t = 0;
	for(std::size_t i = 0; i < 4; i++)
		t += m[i][i];
	return t;
}

X operator*(const M& m, const X& x) {
	X x1;
	for(std::size_t i = 0; i < 4; i++) {
		x1[i] = 0;
		for(std::size_t j = 0; j < 4; j++) {
			x1[i] += m[i][j] * x[j];
		}
	}
	return x1;
}

X operator-(const X& x) {
	X x1;
	for(std::size_t i = 0; i < 4; i++)
		x1[i] = -x[i];
	return x1;
}

X operator+(const X& x1, const X& x2) {
	X x3;
	for(std::size_t i = 0; i < 4; i++)
		x3[i] = x1[i]+x2[i];
	return x3;
}

std::ostream& operator<<(std::ostream& o, const X& x) {
	for(std::size_t i = 0; i < x.size()-1; i++)
		o << x[i] << ", ";
	o << x[x.size()-1];
	return o;
}

ChemicalSystem::ChemicalSystem(
			double _K1, double _K2, double _K3,
			double H, double P, double C,
			double N,
			double S, double SH, double SP, double SC
			) :
	_K1(_K1), _K2(_K2), _K3(_K3),
	H(H), P(P), C(C), N(N), S(S), SH(SH), SP(SP), SC(SC) {
		std::cout << "Init H: " << H/mol << "mol = " << (H/mol)/(V/l) << "mol/L" << std::endl;
		std::cout << "Init P: " << P/mol << "mol = " << (P/mol)/(V/l) << "mol/L = " <<
			((P/entities)*(30.973*u)/mg) / (mineral_weight/kg) << "mg/kg" << std::endl;
		std::cout << "Init C: " << C/mol << "mol = " << (C/mol)/(V/l) << "mol/L = " <<
			((C/entities)*(12.010*u)/mg) / (mineral_weight/kg) << "mg/kg" << std::endl;
		std::cout << std::endl;
		std::cout << "N: " << N << " sites (mol)." << std::endl;
		std::cout << "Init S : " << S/mol << " sites (mol)." << std::endl;
		std::cout << "Init SH: " << SH/mol << " sites (mol)." << std::endl;
		std::cout << "Init SP: " << SP/mol << " sites (mol)." << std::endl;
		std::cout << "Init SC: " << SC/mol << " sites (mol)." << std::endl;
	}

ChemicalSystem ChemicalSystem::equilibrium(
		double K1, double K2, double K3,
		double pH, double solution_P, double solution_C,
		double mineral_N
		) {
	double _K1 = K1/V;
	double _K2 = K2/std::pow(V, 2);
	double _K3 = K3/V;
	double H = V*std::pow(10, -pH)*mol/l;
	double P = ((mineral_weight * solution_P)/(30.973*u)) * entities;
	double C = ((mineral_weight * solution_C)/(12.010*u)) * entities;
	double N = mineral_N * mineral_weight;

	// Assumes that the system starts at equilibrium
	double S = N/(1+_K1*H+_K2*H*P+_K3*C); // S+SH+SP+SC=N
	double SH = _K1 * S * H;
	double SP = _K2 * S * P * H;
	double SC = _K3 * S * C;

	return ChemicalSystem(
			_K1, _K2, _K3,
			H, P, C, N,
			S, SH, SP, SC);
}

ChemicalSystem ChemicalSystem::defaultEquilibrium() {
	return equilibrium(
			std::pow(10, 8.5), // K1
			std::pow(10, 26.3), // K2
			std::pow(10, 0.6), // K3
			7.5, // pH
			1.43e-6 * gram / gram, // 1.43e-6 grams of P in solution by gram of soil
			729.0e-6 * gram / gram, // 729.0e-6 grams of C in solution by gram of soil
			1.45e19 * entities / gram // 1.45e19 sites by gram of soil
			);
}

/*
 *ChemicalSystem ChemicalSystem::soilParameters(
 *            double K1, double K2, double K3,
 *            double pH, double mineral_P, double mineral_C,
 *            double mineral_N
 *            ) {
 *    double _K1 = K1/V;
 *    double _K2 = K2/std::pow(V, 2);
 *    double _K3 = K3/V;
 *    double N = mineral_N * mineral_weight;
 *    double H = V*std::pow(10, -pH)*mol/l;
 *    double P = ((mineral_weight * mineral_P)/(30.973*u)) * entities;
 *    double C = ((mineral_weight * mineral_C)/(12.010*u)) * entities;
 *
 *    ChemicalSystem system(
 *            _K1, _K2, _K3,
 *            H/2, P/2, C/2, N,
 *            S, H/2, P/2, C/2);
 *    system.setPH(pH);
 *    return system;
 *}
 */

/*
 *ChemicalSystem ChemicalSystem::defaultSoil() {
 *    return soilParameters(
 *            std::pow(10, 8.5), // K1
 *            std::pow(10, 26.3), // K2
 *            std::pow(10, 0.6), // K3
 *            7.5, // pH
 *            1.43e-6 * gram / gram, // 1.43e-6 grams of P in solution by gram of soil
 *            729.0e-6 * gram / gram, // 729.0e-6 grams of C in solution by gram of soil
 *            1.45e19 * entities / gram // 1.45e19 sites by gram of soil
 *            );
 *}
 */

double ChemicalSystem::concentration(const double& n) const {
	return (n/V)/(mol/l);
}

double ChemicalSystem::cH() const {
	return concentration(H);
}

double ChemicalSystem::cP() const {
	return concentration(P);
}

double ChemicalSystem::cC() const {
	return concentration(C);
}

double ChemicalSystem::cS() const {
	return S/N;
}

double ChemicalSystem::cSH() const {
	return SH/N;
}

double ChemicalSystem::cSP() const {
	return SP/N;
}

double ChemicalSystem::cSC() const {
	return SC/N;
}

double ChemicalSystem::nH() const {
	return H;
}

double ChemicalSystem::nP() const {
	return P;
}

double ChemicalSystem::nC() const {
	return C;
}

double ChemicalSystem::nS() const {
	return S;
}

double ChemicalSystem::nSH() const {
	return SH;
}

double ChemicalSystem::nSP() const {
	return SP;
}

double ChemicalSystem::nSC() const {
	return SC;
}

X ChemicalSystem::reactionQuotient() const {
	return {{
		0,
		SH / (S*H),
		SP / (S*H*P),
		SC / (S*C)
	}};
}

void ChemicalSystem::distanceToEquilibrium() const {
	X reaction_quotient = reactionQuotient();
	std::cout << "pH: " << -std::log10((H/V)/(mol/l)) << std::endl;
	std::cout << "K1: " << reaction_quotient[1]/_K1 << std::endl;
	std::cout << "K2: " << reaction_quotient[2]/_K2 << std::endl;
	std::cout << "K3: " << reaction_quotient[3]/_K3 << std::endl;
}

void ChemicalSystem::incrementP(double dP) {
	P += dP;
}

void ChemicalSystem::setPH(double pH) {
	double alpha = V*std::pow(10, -pH)*mol/l - H;
	std::cout << "alpha: " << alpha << std::endl;

	F f(alpha, *this);
	X dx = Newton<X, M>(
			{{0, 0, 0, 0}},
			[&f] (const X& x) {return f.f(x);},
			[&f] (const X& x) {return f.df(x);}
			).solve_iter(10);
	S = S - dH(dx) - dP(dx) - dC(dx);
	H = H + deltaH(dx) - dH(dx) - dP(dx);
	P = P  - dP(dx);
	C = C - dC(dx);
	SH = SH + dH(dx);
	std::cout << "dP(dx): " << dP(dx) << std::endl;
	SP = SP + dP(dx);
	SC = SC + dC(dx);
	std::cout << "alpha end: " << V*std::pow(10, -pH)*mol/l - H << std::endl;
}

X ChemicalSystem::F::f(const X& x) {
	return {{
		deltaH(x) - (dH(x)+dP(x)) - alpha,
		(problem.SH + dH(x)) /
			((problem.S-dH(x)-dP(x)-dC(x))*(problem.H+deltaH(x)-dH(x)-dP(x)))
			/ problem._K1 - 1,
		(problem.SP + dP(x)) /
			((problem.S-dH(x)-dP(x)-dC(x))*(problem.H+deltaH(x)-dH(x)-dP(x))*(problem.P-dP(x)))
			/ problem._K2 - 1,
		(problem.SC + dC(x)) /
			((problem.S-dH(x)-dP(x)-dC(x))*(problem.C-dC(x)))
			/ problem._K3 - 1
	}};
}

M ChemicalSystem::F::df(const X& x) {
	// utils
	double v;
	double v2;

	double DS = problem.S-dH(x)-dP(x)-dC(x);
	double DH = problem.H+deltaH(x)-dH(x)-dP(x);
	double DP = problem.P-dP(x);
	double DC = problem.C-dC(x);
	std::cout << "  DP: " << DP << ", DH: " << DH << ", DS: " << DS << ", DC: " << DC << std::endl;

	// f1
	double df1_dx1 = 1;
	double df1_dx2 = -1;
	double df1_dx3 = -1;
	double df1_dx4 = 0;

	std::cout << "df1" << std::endl;
	std::cout << "  df1_dx2: " << df1_dx2 << std::endl;

	// f2
	v = DS*DH;
	v2 = std::pow(v, 2);
	double df2_dx1 = 1/problem._K1*(-DS*(problem.SH + dH(x))/v2);
	double df2_dx2 = 1/problem._K1*((DS*DH + (DH+DS)*(problem.SH+dH(x)))/v2);
	double df2_dx3 = 1/problem._K1*(((DH+DS)*(problem.SH+dH(x)))/v2);
	double df2_dx4 = 1/problem._K1*(DH*(problem.SH+dH(x))/v2);
	std::cout << "df3" << std::endl;
	std::cout << "  v2: " << v2 << std::endl;
	std::cout << "  df2_dx2: " << df2_dx2 << std::endl;

	// f3
	v = DS*DH*DP;
	v2 = std::pow(v, 2);
	double df3_dx1 = 1/problem._K2*(-DP*DS*(problem.SP+dP(x))/v2);
	double df3_dx2 = 1/problem._K2*(DP*(DH+DS)*(problem.SP+dP(x))/v2);
	double df3_dx3 = 1/problem._K2*((v+(DS*DH + DP*DH + DS*DP)*(problem.SP+dP(x)))/v2);
	double df3_dx4 = 1/problem._K2*(DH*DP*(problem.SP+dP(x))/v2);
	std::cout << "df3" << std::endl;
	std::cout << "  v2: " << v2 << std::endl;
	std::cout << "  df3_dx2: " << df3_dx2 << std::endl;

	// f4
	v = DS*DC;
	v2 = std::pow(v, 2);
	double df4_dx1 = 0;
	double df4_dx2 = 1/problem._K3*(DC * (problem.SC+dC(x))/v2);
	double df4_dx3 = df4_dx2;
	double df4_dx4 = 1/problem._K3*((DS*DC + (DC+DS)*(problem.SC+dC(x)))/v2);
	std::cout << "df4" << std::endl;
	std::cout << "  v2: " << v2 << std::endl;
	std::cout << "  df4_dx2: " << df2_dx2 << std::endl;
	
	return {{
		{{df1_dx1, df1_dx2, df1_dx3, df1_dx4}},
		{{df2_dx1, df2_dx2, df2_dx3, df2_dx4}},
		{{df3_dx1, df3_dx2, df3_dx3, df3_dx4}},
		{{df4_dx1, df4_dx2, df4_dx3, df4_dx4}}
	}};
}

