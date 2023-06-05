#include <array>
#include <cmath>
#include <iostream>
#include "linear.h"
#include "units.h"

namespace mineral {
	class ChemicalSystem;

	namespace ph_solver {
		typedef mineral::X<double, 4> X;
		typedef mineral::M<double, 4> M;

		double deltaH(const X& x);
		double dH(const X& x);
		double dP(const X& x);
		double dC(const X& x);

		class RawPHSolver {
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

			public:
			static X solve(const ChemicalSystem& system, double pH);

		};

		class LogPHSolver {
			class F {
				double pH;
				const ChemicalSystem& system;

				public:
				F(double pH, const ChemicalSystem& system)
					: pH(pH), system(system) {
					}
				X f(const X& x);
				M df(const X& x);
			};

			public:
			static X solve(const ChemicalSystem& system, double pH);
		};
	}

	namespace equilibrium_solver {
		typedef mineral::X<double, 3> X;
		typedef mineral::M<double, 3> M;

		double dH(const X& x);
		double dP(const X& x);
		double dC(const X& x);

		class Solver {
			class F {
				const ChemicalSystem& system;

				public:
				F(const ChemicalSystem& system) : system(system) {
					}
				X f(const X& x);
				M df(const X& x);
			};

			public:
			static X solve(const ChemicalSystem& system);
		};
	}

	class ChemicalSystem {
		public:
		static const double bulk_density;
		static const double V;
		static const double mineral_weight;

		private:
		double _K1;
		double _K2;
		double _K3;
		double _initPH;

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
				double K1, double K2, double K3, double initPH,
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
		static ChemicalSystem Devau2011Control();

		//static ChemicalSystem soilParameters(
		//double K1, double K2, double K3,
		//double pH, double soil_P, double soil_C,
		//double mineral_N
		//);
		//static ChemicalSystem defaultSoil();


		void incrementP(double P);
		void setPH(double pH);

		void solveEquilibrium();

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
		double initPH() const;

		double pH() const;
	};
}
