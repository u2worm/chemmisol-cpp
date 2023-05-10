#include <array>
#include <cmath>
#include <iostream>
#include "linear.h"

namespace mineral {
	class ChemicalSystem;

	namespace pHSolver {
		typedef mineral::X<double, 4> X;
		typedef mineral::M<double, 4> M;

		double deltaH(const X& x);
		double dH(const X& x);
		double dP(const X& x);
		double dC(const X& x);

		class F {
			double alpha;
			const ChemicalSystem& system;

			public:
			F(double alpha, const ChemicalSystem& system)
				: alpha(alpha), system(system) {
				}
			X f(const X& x);
			M df(const X& x);
		};

		X solve(const ChemicalSystem& system, double pH);
	}

	static const double m = 1;
	static const double cm = 1e-2 * m;
	static const double m3 = m*m*m;
	static const double cm3 = cm*cm*cm;
	static const double l = 1e-3 * m3;

	static const double gram = 1;
	static const double kg = 1e3 * gram;
	static const double mg = 1e-3 * gram;

	static const double NA = 6.02214076e23;
	static const double mol = 1;
	static const double entities = 1/NA;

	static const double u = 1.66053906660e-27*kg;

	class ChemicalSystem {
		public:
		static const double bulk_density;
		static const double V;
		static const double mineral_weight;

		private:
		double _K1;
		double _K2;
		double _K3;

		double H;
		double P;
		double C;

		double N;
		double S;
		double SH;
		double SP;
		double SC;

		double concentration(const double& n) const;

		private:
		ChemicalSystem(
				double K1, double K2, double K3,
				double H, double P, double C,
				double N,
				double S, double SH, double SP, double SC
				);
		public:
		
		static ChemicalSystem equilibrium(
				double K1, double K2, double K3,
				double pH, double solution_P, double solution_C,
				double mineral_N
				);
		static ChemicalSystem defaultEquilibrium();

		//static ChemicalSystem soilParameters(
		//double K1, double K2, double K3,
		//double pH, double soil_P, double soil_C,
		//double mineral_N
		//);
		//static ChemicalSystem defaultSoil();


		void incrementP(double P);
		void setPH(double pH);

		//X reactionQuotient() const;
		void distanceToEquilibrium() const;

		double cH() const;
		double cP() const;
		double cC() const;
		double cS() const;
		double cSH() const;
		double cSP() const;
		double cSC() const;

		double nH() const;
		double nP() const;
		double nC() const;
		double nS() const;
		double nSH() const;
		double nSP() const;
		double nSC() const;

		double K1() const;
		double K2() const;
		double K3() const;
	};
}
