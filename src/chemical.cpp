#include "chemical.h"
#include "newton.h"

#include <iostream>
#include <limits>

namespace mineral {
	const double ChemicalSystem::bulk_density = 1.17*gram/cm3;
	const double ChemicalSystem::V = std::pow(1*cm, 3);
	const double ChemicalSystem::mineral_weight  = V*bulk_density;

	namespace ph_solver {
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

		X RawPHSolver::F::f(const X& x) {
			return X({
					deltaH(x) - (dH(x)+dP(x)) - alpha,
					(system.nSH() + dH(x)) /
					((system.nS()-dH(x)-dP(x)-dC(x))*(system.nH()+deltaH(x)-dH(x)-dP(x)))
					/ system.K1() - 1,
					(system.nSP() + dP(x)) /
					((system.nS()-dH(x)-dP(x)-dC(x))*(system.nH()+deltaH(x)-dH(x)-dP(x))*(system.nP()-dP(x)))
					/ system.K2() - 1,
					(system.nSC() + dC(x)) /
					((system.nS()-dH(x)-dP(x)-dC(x))*(system.nC()-dC(x)))
					/ system.K3() - 1
					});
		}

		M RawPHSolver::F::df(const X& x) {
			// utils
			double v;
			double v2;

			double DH = system.nH()+deltaH(x)-dH(x)-dP(x);
			double DP = system.nP()-dP(x);
			double DC = system.nC()-dC(x);
			double DS = system.nS()-dH(x)-dP(x)-dC(x);
			double DSH = system.nSH() + dH(x);
			double DSP = system.nSP() + dP(x);
			double DSC = system.nSC() + dC(x);
			std::cout << "  DP: " << DP << ", DH: " << DH << ", DC: " << DC << std::endl;
			std::cout << "  DS: " << DS << ", DSH: " << DSH << ", DSP: " << DSP << ", DSC: " << DSC << std::endl;

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
			double df2_dx1 = 1/system.K1()*(-DS*DSH/v2);
			double df2_dx2 = 1/system.K1()*((DS*DH + (DH+DS)*DSH)/v2);
			double df2_dx3 = 1/system.K1()*(((DH+DS)*DSH)/v2);
			double df2_dx4 = 1/system.K1()*(DH*DSH/v2);
			std::cout << "df3" << std::endl;
			std::cout << "  v2: " << v2 << std::endl;
			std::cout << "  df2_dx2: " << df2_dx2 << std::endl;

			// f3
			v = DS*DH*DP;
			v2 = std::pow(v, 2);
			double df3_dx1 = 1/system.K2()*(-DP*DS*DSP/v2);
			double df3_dx2 = 1/system.K2()*(DP*(DH+DS)*DSP/v2);
			double df3_dx3 = 1/system.K2()*((v+(DS*DH + DP*DH + DS*DP)*DSP)/v2);
			double df3_dx4 = 1/system.K2()*(DH*DP*DSP/v2);
			std::cout << "df3" << std::endl;
			std::cout << "  v2: " << v2 << std::endl;
			std::cout << "  df3_dx2: " << df3_dx2 << std::endl;

			// f4
			v = DS*DC;
			v2 = std::pow(v, 2);
			double df4_dx1 = 0;
			double df4_dx2 = 1/system.K3()*(DC * DSC/v2);
			double df4_dx3 = df4_dx2;
			double df4_dx4 = 1/system.K3()*((DS*DC + (DC+DS)*DSC)/v2);
			std::cout << "df4" << std::endl;
			std::cout << "  v2: " << v2 << std::endl;
			std::cout << "  df4_dx2: " << df2_dx2 << std::endl;

			return M({
					{df1_dx1, df1_dx2, df1_dx3, df1_dx4},
					{df2_dx1, df2_dx2, df2_dx3, df2_dx4},
					{df3_dx1, df3_dx2, df3_dx3, df3_dx4},
					{df4_dx1, df4_dx2, df4_dx3, df4_dx4}
					});
		}

		X RawPHSolver::solve(const ChemicalSystem& system, double pH) {
			double alpha = ChemicalSystem::V*std::pow(10, -pH)*mol/l - system.nH();
			std::cout << "alpha: " << alpha << std::endl;

			F f(alpha, system);
			return Newton<X, M>(
					{{0, 0, 0, 0}},
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(20);
		}

		X LogPHSolver::F::f(const X& x) {
			double DH = system.nH()+deltaH(x)-dH(x)-dP(x);
			double DP = system.nP()-dP(x);
			double DC = system.nC()-dC(x);
			double DS = system.nS()-dH(x)-dP(x)-dC(x);
			double DSH = system.nSH() + dH(x);
			double DSP = system.nSP() + dP(x);
			double DSC = system.nSC() + dC(x);
			std::cout << "  DP: " << DP << ", DH: " << DH << ", DC: " << DC << std::endl;
			std::cout << "  DS: " << DS << ", DSH: " << DSH << ", DSP: " << DSP << ", DSC: " << DSC << std::endl;
			return X({
					log(DH) - log(ChemicalSystem::V) + pH - log(mol/l),
					log(DSH)-log(DS)-log(DH)-log(system.K1()),
					log(DSP)-log(DS)-log(DH)-log(DP)-log(system.K2()),
					log(DSC)-log(DS)-log(DC)-log(system.K3())
					});
		}

		M LogPHSolver::F::df(const X& x) {
			// utils
			//double v;
			//double v2;

			double DH = system.nH()+deltaH(x)-dH(x)-dP(x);
			double DP = system.nP()-dP(x);
			double DC = system.nC()-dC(x);
			double DS = system.nS()-dH(x)-dP(x)-dC(x);
			double DSH = system.nSH() + dH(x);
			double DSP = system.nSP() + dP(x);
			double DSC = system.nSC() + dC(x);

			// f1
			double df1_dx1 = (1/DH)/ln10;
			double df1_dx2 = -(1/DH)/ln10;
			double df1_dx3 = -(1/DH)/ln10;
			double df1_dx4 = 0;

			//std::cout << "df1" << std::endl;
			//std::cout << "  df1_dx2: " << df1_dx2 << std::endl;

			// f2
			//v = DS*DH;
			//v2 = std::pow(v, 2);
			double df2_dx1 = -(1/DH)/ln10;
			double df2_dx2 = (1/DSH + 1/DS + 1/DH)/ln10;
			double df2_dx3 = (1/DS + 1/DH)/ln10;
			double df2_dx4 = (1/DS)/ln10;
			//std::cout << "df3" << std::endl;
			//std::cout << "  v2: " << v2 << std::endl;
			//std::cout << "  df2_dx2: " << df2_dx2 << std::endl;

			// f3
			//v = DS*DH*DP;
			//v2 = std::pow(v, 2);
			double df3_dx1 = df2_dx1;
			double df3_dx2 = (1/DS + 1/DH)/ln10;
			double df3_dx3 = (1/DSP + 1/DS + 1/DH + 1/DP)/ln10;
			double df3_dx4 = (1/DS)/ln10;
			//std::cout << "df3" << std::endl;
			//std::cout << "  v2: " << v2 << std::endl;
			//std::cout << "  df3_dx2: " << df3_dx2 << std::endl;

			// f4
			//v = DS*DC;
			//v2 = std::pow(v, 2);
			double df4_dx1 = 0;
			double df4_dx2 = (1/DS)/ln10;
			double df4_dx3 = df4_dx2;
			double df4_dx4 = (1/DSC + 1/DS + 1/DC)/ln10;
			//std::cout << "df4" << std::endl;
			//std::cout << "  v2: " << v2 << std::endl;
			//std::cout << "  df4_dx2: " << df2_dx2 << std::endl;

			return M({
					{df1_dx1, df1_dx2, df1_dx3, df1_dx4},
					{df2_dx1, df2_dx2, df2_dx3, df2_dx4},
					{df3_dx1, df3_dx2, df3_dx3, df3_dx4},
					{df4_dx1, df4_dx2, df4_dx3, df4_dx4}
					});
		}

		X LogPHSolver::solve(const ChemicalSystem& system, double pH) {
			F f(pH, system);
			return Newton<X, M>(
					{{0, 0, 0, 0}},
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(20);
		}
	}

	namespace equilibrium_solver {
		double dH(const X& x) {
			return x[0];
		}
		double dP(const X& x) {
			return x[1];
		}
		double dC(const X& x) {
			return x[2];
		}

		X Solver::F::f(const X &x) {
			double DH = system.nH()-dH(x)-dP(x);
			double DP = system.nP()-dP(x);
			double DC = system.nC()-dC(x);
			double DS = system.nS()-dH(x)-dP(x)-dC(x);
			double DSH = system.nSH() + dH(x);
			double DSP = system.nSP() + dP(x);
			double DSC = system.nSC() + dC(x);
			std::cout << "  DP: " << DP << ", DH: " << DH << ", DC: " << DC << std::endl;
			std::cout << "  DS: " << DS << ", DSH: " << DSH << ", DSP: " << DSP << ", DSC: " << DSC << std::endl;
			return X({
					log(DSH)-log(DS)-log(DH)-log(system.K1()),
					log(DSP)-log(DS)-log(DH)-log(DP)-log(system.K2()),
					log(DSC)-log(DS)-log(DC)-log(system.K3())
					});
		}

		M Solver::F::df(const X& x) {
			double DH = system.nH()-dH(x)-dP(x);
			double DP = system.nP()-dP(x);
			double DC = system.nC()-dC(x);
			double DS = system.nS()-dH(x)-dP(x)-dC(x);
			double DSH = system.nSH() + dH(x);
			double DSP = system.nSP() + dP(x);
			double DSC = system.nSC() + dC(x);

			// f1
			double df1_dx1 = (1/DSH + 1/DS + 1/DH)/ln10;
			double df1_dx2 = (1/DS + 1/DH)/ln10;
			double df1_dx3 = (1/DS)/ln10;

			// f2
			double df2_dx1 = (1/DS + 1/DH)/ln10;
			double df2_dx2 = (1/DSP + 1/DS + 1/DH + 1/DP)/ln10;
			double df2_dx3 = (1/DS)/ln10;

			// f3
			double df3_dx1 = (1/DS)/ln10;
			double df3_dx2 = df3_dx1;
			double df3_dx3 = (1/DSC + 1/DS + 1/DC)/ln10;

			return M({
					{df1_dx1, df1_dx2, df1_dx3},
					{df2_dx1, df2_dx2, df2_dx3},
					{df3_dx1, df3_dx2, df3_dx3}
					});
		}

		X Solver::solve(const ChemicalSystem& system) {
			F f(system);
			return Newton<X, M>(
					{{0, 0, 0, 0}},
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(20);
		}

	}


	ChemicalSystem::ChemicalSystem(
			double _K1, double _K2, double _K3, double initPH,
			double H, double P, double C,
			double N,
			double S, double SH, double SP, double SC
			) :
		_K1(_K1), _K2(_K2), _K3(_K3), _initPH(initPH),
		H(H), P(P), C(C), N(N), S(S), SH(SH), SP(SP), SC(SC) {
			std::cout << "[Init] Mineral weight: " << mineral_weight/kg << "kg" << std::endl;
			std::cout << "[Init] H: " << H/mol << "mol = " << (H/mol)/(V/l) << "mol/L" << std::endl;
			std::cout << "[Init] P: " << P/mol << "mol = " << (P/mol)/(V/l) << "mol/L = " <<
				((P/entities)*(30.973*u)/mg) / (mineral_weight/kg) << "mg/kg" << std::endl;
			std::cout << "[Init] C: " << C/mol << "mol = " << (C/mol)/(V/l) << "mol/L = " <<
				((C/entities)*(12.010*u)/mg) / (mineral_weight/kg) << "mg/kg" << std::endl;
			std::cout << "[Init] pH: " << this->initPH() << std::endl;
			std::cout << std::endl;
			std::cout << "[Init] N: " << N << " sites (mol)." << std::endl;
			std::cout << "[Init] S : " << S/mol << " sites (mol)." << std::endl;
			std::cout << "[Init] SH: " << SH/mol << " sites (mol)." << std::endl;
			std::cout << "[Init] SP: " << SP/mol << " sites (mol)." << std::endl;
			std::cout << "[Init] SC: " << SC/mol << " sites (mol)." << std::endl;
		}

	ChemicalSystem ChemicalSystem::equilibrium(
			double K1, double K2, double K3, double pH,
			double solution_P, double solution_C,
			double mineral_N
			) {
		double _K1 = K1/V;
		double _K2 = K2/std::pow(V, 2);
		double _K3 = K3/V;
		double H = V*std::pow(10, -pH)*mol/l;
		double P = ((mineral_weight * solution_P)/(30.973*u)) * entities;
		double C = ((mineral_weight * solution_C)/(12.010*u)) * entities;
		double N = mineral_N;

		// Assumes that the system starts at equilibrium
		double S = N/(1+_K1*H+_K2*H*P+_K3*C); // S+SH+SP+SC=N
		double SH = _K1 * S * H;
		double SP = _K2 * S * P * H;
		double SC = _K3 * S * C;

		return ChemicalSystem(
				_K1, _K2, _K3, pH,
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
				1.45e19 * entities / gram * mineral_weight // 1.45e19 sites by gram of soil
				);
	}

	ChemicalSystem ChemicalSystem::Devau2011Control() {
		return equilibrium(
				// Kaolinite
				std::pow(10, 4.36),
				std::pow(10, 23),
				std::pow(10, 1),
				6.5,
				(279-69)*mg/kg,
				9.80*gram/kg,
				(6.15*entities/nm2)*(105*m2/gram)*(20.12*gram/l)*V
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

	double ChemicalSystem::K1() const {
		return _K1;
	}

	double ChemicalSystem::K2() const {
		return _K2;
	}

	double ChemicalSystem::K3() const {
		return _K3;
	}

	double ChemicalSystem::initPH() const {
		return _initPH;
	}

	double ChemicalSystem::pH() const {
		return -log((nH()/V)/(mol/l));
	}

	/*
	 *X ChemicalSystem::reactionQuotient() const {
	 *    return {{
	 *        0,
	 *        SH / (S*H),
	 *        SP / (S*H*P),
	 *        SC / (S*C)
	 *    }};
	 *}
	 *
	 *void ChemicalSystem::distanceToEquilibrium() const {
	 *    X reaction_quotient = reactionQuotient();
	 *    std::cout << "pH: " << -std::log10((H/V)/(mol/l)) << std::endl;
	 *    std::cout << "K1: " << reaction_quotient[1]/_K1 << std::endl;
	 *    std::cout << "K2: " << reaction_quotient[2]/_K2 << std::endl;
	 *    std::cout << "K3: " << reaction_quotient[3]/_K3 << std::endl;
	 *}
	 */

	void ChemicalSystem::incrementP(double dP) {
		P += dP;
	}

	void ChemicalSystem::setPH(double pH) {
		std::cout << "[PH SOLVER] start" << std::endl;
		using namespace ph_solver;
		ph_solver::X dx = LogPHSolver::solve(*this, pH);
		S = S - dH(dx) - dP(dx) - dC(dx);
		H = H + deltaH(dx) - dH(dx) - dP(dx);
		P = P - dP(dx);
		C = C - dC(dx);
		SH = SH + dH(dx);
		std::cout << "deltaH(dx): " << deltaH(dx) << std::endl;
		std::cout << "dH(dx): " << dH(dx) << std::endl;
		std::cout << "dP(dx): " << dP(dx) << std::endl;
		std::cout << "dC(dx): " << dC(dx) << std::endl;
		SP = SP + dP(dx);
		SC = SC + dC(dx);
		std::cout << "pH end: " << -log((H/V)/(mol/l)) << std::endl;
		std::cout << "Q1    : " << SH/(S*H)/_K1 << std::endl;
		std::cout << "Q2    : " << SP/(S*P*H)/_K2 << std::endl;
		std::cout << "Q3    : " << SC/(S*C)/_K3 << std::endl;

		std::cout << "[PH SOLVER] end" << std::endl;
	}

	void ChemicalSystem::solveEquilibrium() {
		std::cout << "[EQ SOLVER] start" << std::endl;
		using namespace equilibrium_solver;
		equilibrium_solver::X dx = Solver::solve(*this);
		S = S - dH(dx) - dP(dx) - dC(dx);
		H = H - dH(dx) - dP(dx);
		P = P - dP(dx);
		C = C - dC(dx);
		SH = SH + dH(dx);
		std::cout << "dH(dx): " << dH(dx) << std::endl;
		std::cout << "dP(dx): " << dP(dx) << std::endl;
		std::cout << "dC(dx): " << dC(dx) << std::endl;
		SP = SP + dP(dx);
		SC = SC + dC(dx);
		std::cout << "pH end: " << -log((H/V)/(mol/l)) << std::endl;
		std::cout << "Q1    : " << SH/(S*H)/_K1 << std::endl;
		std::cout << "Q2    : " << SP/(S*P*H)/_K2 << std::endl;
		std::cout << "Q3    : " << SC/(S*C)/_K3 << std::endl;

		std::cout << "[EQ SOLVER] end" << std::endl;
	}
}

