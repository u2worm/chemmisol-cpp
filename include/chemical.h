#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <memory>
#include "linear.h"
#include "units.h"

namespace mineral {
	class ChemicalSystem;

	namespace solver {
		typedef std::vector<double> X;
		typedef std::vector<std::vector<double>> M;

		class Solver {
			public:
				class F {
					const ChemicalSystem& system;
					private:
					std::vector<double> init_concentrations;

					public:
					X concentrations(const X& extent) const;
					F(const ChemicalSystem& system);
					X f(const X& extent) const;
					M df(const X& extent) const;
				};

				static X solve(const ChemicalSystem& system);
		};
	}

	class Component {
		public:
			static const double V;
		private:
			std::string name;
			std::size_t id;

		public:
			Component(const std::string& name, std::size_t id)
				: name(name), id(id) {
				}

			const std::string& getName() const {
				return name;
			}

			std::size_t getId() const {
				return id;
			}

			virtual void incrementConcentration(double extent) = 0;
			virtual double concentration(double current_concentration, double extent) const = 0;
			virtual double concentration() const = 0;

			double activity() const {
				return activity(concentration());
			}

			virtual double activity(double current_concentration) const = 0;

			/**
			 * Computes the activity that results from an increment of
			 * "extent" of the quantity of the current component.
			 *
			 * @param current_concentrations current concentrations of all
			 * component. The concentration of the current component can be
			 * retrieved with current_concentrations[getId()].
			 * @param extent increment of the quantity of the current component
			 * (mole)
			 */
			virtual double activity(
					const std::vector<double>& current_concentrations,
					double extent) const = 0;

			/**
			 * Computes dc_k/dx_j where c_k denotes the concentration of the
			 * current component, and x_j denotes the extent of reaction j.
			 *
			 * @param current_concentrations Current concentrations of all
			 * components. The concentration of the current component can be
			 * retrieved with current_concentrations[getId()].
			 * @param d_coef Stoichiometric coefficient of the current component
			 * in the reaction k. Might be 0.
			 */
			virtual double Dactivity(
					const std::vector<double>& current_concentrations,
					double d_coef
					) const = 0;
			
			virtual ~Component() {
			}
	};

	class Solvent : public Component {
		public:
			Solvent(const std::string& name, std::size_t id)
				: Component(name, id) {
				}

			void incrementConcentration(double) override {
			}

			double concentration(double, double) const override {
				return 0.0;
			}

			double concentration() const override {
				return 1.0;
			}

			double activity(double) const override {
				return 1.0;
			}

			double activity(const std::vector<double>&,  double) const override {
				return 1.0;
			}

			double Dactivity(
					const std::vector<double>&,
					double
					) const override {
				return 0.0;
			}
	};

	class AqueousComponent : public Component {
		private:
			double C;
		
		public:
			AqueousComponent(
					const std::string& name, std::size_t id,
					double C)
				: Component(name, id), C(C) {
				}

			void incrementConcentration(double extent) override {
				C+=extent/V;
			}

			double concentration(double current_concentration, double extent) const override{
				return current_concentration + extent/V;
			}

			double concentration() const override {
				return C;
			}

			double activity(double current_concentration) const override {
				return current_concentration/(1*mol/l);
			}

			double activity(
					const std::vector<double>& current_concentrations,
					double extent) const override {
				return activity(current_concentrations[getId()] + (extent/V));
			}

			double Dactivity(
					const std::vector<double>&,
					double d_coef
					) const override {
				// Note: this is 0 if d_coef=0, i.e. if the reaction k does not
				// use this component.
				return (d_coef/V)/(1*mol/l);
			}
	};

	class MineralComponent : public Component {
		private:
			double fraction; // S/N
			double N;

		public:
			MineralComponent(
					const std::string& name, std::size_t id,
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					double fraction)
				: Component(name, id), fraction(fraction),
				N(V*solid_concentration*specific_surface_area*site_concentration) {
				}
			MineralComponent(
					const std::string& name, std::size_t id,
					double fraction)
				: Component(name, id), fraction(fraction) {
				}

			void incrementConcentration(double extent) override {
				fraction+=extent/N;
			}

			double concentration(double current_concentration, double extent) const override{
				return current_concentration + extent/N;
			}
			double concentration() const override {
				return fraction;
			}

			double activity(double current_concentration) const override {
				return current_concentration;
			}

			double activity(
					const std::vector<double>& current_concentrations,
					double extent) const override {
				return activity(current_concentrations[getId()]+(extent/N));
			}

			double Dactivity(
					const std::vector<double>&,
					double d_coef
					) const override {
				// Note: this is 0 if d_coef=0, i.e. if the reaction k does not
				// use this component.
				return d_coef/N;
			}

	};

	class Reaction {
		private:
			std::string name;
			std::size_t id;
			double log_K;
			std::unordered_map<std::string, double> reactives;

		public:
			Reaction(
					const std::string& name, std::size_t id, double log_K,
					const std::unordered_map<std::string, double>& reactives
					)
				: name(name), id(id), log_K(log_K), reactives(reactives) {
				}

			const std::string& getName() const {
				return name;
			}

			std::size_t getId() const {
				return id;
			}

			double getLogK() const {
				return log_K;
			}

			const std::unordered_map<std::string, double> getReactives() const {
				return reactives;
			}
	};

	class ChemicalSystem {
		private:
			std::size_t component_index = 0;
			std::size_t reaction_index = 0;

			std::vector<std::unique_ptr<Component>> components;
			std::unordered_map<std::string, const Component*> components_by_name;
			std::vector<std::unique_ptr<Reaction>> reactions;
			std::unordered_map<std::string, const Reaction*> reactions_by_name;
			std::vector<std::vector<double>> reaction_matrix;

			void addComponent(Component* component);
			void addReaction(Reaction* reaction);

			std::size_t max_iteration = 20;
		public:
			void addReaction(
					std::string name,
					double log_K,
					std::unordered_map<std::string, double> reactives
					);
			void addComponent(
					std::string name,
					double concentration
					);
			void initPH(double pH);
			void fixPH(double pH);

			const Reaction& getReaction(const std::string& name) const;
			const Component& getComponent(const std::string& name) const;
			const Component& getComponent(const std::size_t& id) const;

			const std::vector<std::unique_ptr<Component>>& getComponents() const {
				return components;
			}

			const std::vector<std::unique_ptr<Reaction>>& getReactions() const {
				return reactions;
			}

			const std::vector<std::vector<double>>& getReactionMatrix() const {
				return reaction_matrix;
			}

			void initReactionMatrix();
			void solveEquilibrium();

			std::size_t getMaxIteration() const {
				return max_iteration;
			}

			void setMaxIteration(const std::size_t& max_iteration) {
				this->max_iteration = max_iteration;
			}
	};

/*
 *    class Soil {
 *        public:
 *        static const double bulk_density;
 *        static const double V;
 *        static const double mineral_weight;
 *
 *        private:
 *
 *        double _K1;
 *        double _K2;
 *        double _K3;
 *        double _initPH;
 *
 *        double H;
 *        double P;
 *        double C;
 *
 *        double N;
 *        double S;
 *        double SH;
 *        double SP;
 *        double SC;
 *
 *        double concentration(const double& n) const;
 *
 *        private:
 *        ChemicalSystem(
 *                double K1, double K2, double K3, double initPH,
 *                double H, double P, double C,
 *                double N,
 *                double S, double SH, double SP, double SC
 *                );
 *
 *        public:
 *        static ChemicalSystem equilibrium(
 *                double K1, double K2, double K3,
 *                double pH, double solution_P, double solution_C,
 *                double mineral_N
 *                );
 *        static ChemicalSystem defaultEquilibrium();
 *        static ChemicalSystem Devau2011Control();
 *
 *        //static ChemicalSystem soilParameters(
 *        //double K1, double K2, double K3,
 *        //double pH, double soil_P, double soil_C,
 *        //double mineral_N
 *        //);
 *        //static ChemicalSystem defaultSoil();
 *
 *
 *        void incrementP(double P);
 *        void setPH(double pH);
 *
 *        //X reactionQuotient() const;
 *        void distanceToEquilibrium() const;
 *
 *        double cH() const;
 *        double cP() const;
 *        double cC() const;
 *        double cS() const;
 *        double cSH() const;
 *        double cSP() const;
 *        double cSC() const;
 *
 *        double nH() const;
 *        double nP() const;
 *        double nC() const;
 *        double nS() const;
 *        double nSH() const;
 *        double nSP() const;
 *        double nSC() const;
 *
 *        double K1() const;
 *        double K2() const;
 *        double K3() const;
 *        double initPH() const;
 *
 *        double pH() const;
 *    };
 */
}
