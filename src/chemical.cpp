#include "chemical.h"
#include "newton.h"

#include <iostream>
#include <limits>
#include <stdexcept>

namespace mineral {
	/*
	 *const double ChemicalSystem::bulk_density = 1.17*gram/cm3;
	 *const double ChemicalSystem::V = std::pow(1*cm, 3);
	 *const double ChemicalSystem::mineral_weight  = V*bulk_density;
	 */

	const double Component::V = 1*l;

	void ChemicalSystem::addReaction(Reaction* reaction) {
		reactions.emplace_back(reaction);
		reactions_by_name.emplace(reaction->getName(), reaction);
	}

	void ChemicalSystem::addReaction(
			std::string name,
			double K,
			std::unordered_map<std::string, double> reactives
			) {
		addReaction(new Reaction(
				name, reaction_index++, K, reactives
				));
	}

	void ChemicalSystem::addComponent(Component* component) {
		components.emplace_back(component);
		components_by_name.emplace(component->getName(), component);
	}
	void ChemicalSystem::addComponent(std::string name, double concentration) {
		addComponent(
				new AqueousComponent(name, component_index++, concentration)
				);
	}

	void ChemicalSystem::initPH(double pH) {
		addComponent("H+", std::pow(10, -pH)*mol/l);
		addComponent("OH-", std::pow(10, -pH)*mol/l);
	}
	void ChemicalSystem::fixPH(double pH) {
		initPH(pH);
		// A fake reaction that produces (or consumes) as much H+ as needed to
		// fix the pH to the provided value
		addReaction("pH", -pH, {{"H+", -1}});
	}

	const Reaction& ChemicalSystem::getReaction(const std::string& name) const {
		return *reactions_by_name.find(name)->second;
	}

	const Component& ChemicalSystem::getComponent(const std::string& name) const {
		return *components_by_name.find(name)->second;
	}

	const Component& ChemicalSystem::getComponent(const std::size_t& id) const {
		return *components[id];
	}

	void ChemicalSystem::initReactionMatrix() {
		for(auto& reaction : reactions) {
			for(auto& component : reaction->getReagents()) {
				if(components_by_name.find(component.first) == components_by_name.end()) {
					if(component.first == "H2O") {
						addComponent(new Solvent("H2O", component_index++));
					} else {
						addComponent(component.first, 0);
					}
				}
			}
		}

		reaction_matrix.resize(reactions.size());
		for(const auto& reaction : reactions) {
			reaction_matrix[reaction->getIndex()].resize(components.size());
			for(const auto& reactive : reaction->getReagents()) {
				reaction_matrix[reaction->getIndex()][getComponent(reactive.first).getIndex()]
					= reactive.second;
			}
		}
	}

	void ChemicalSystem::solveEquilibrium() {
		initReactionMatrix();
		using namespace solver;
		solver::X dx = solve(*this);
		for(auto& reaction : reactions) {
			for(auto& component : reaction->getReagents()) {
				components[getComponent(component.first).getIndex()]->incrementConcentration(
						component.second * dx[reaction->getIndex()]
						);
			}
		}
	}

	double ChemicalSystem::reactionQuotient(const std::string& name) const {
		double quotient = 1;
		auto& reaction = getReaction(name);
		for(auto& reactive : reaction.getReagents()) {
			quotient *= std::pow(getComponent(reactive.first).activity(), reactive.second);
		}
		// 1/quotient due to the following conventions:
		// - products are given with **negative** stoichiometric coefficients
		// - products from the numerous of the reaction quotient
		return 1/quotient;
	}

	namespace solver {
		F::F(const ChemicalSystem& system)
			: system(system), init_concentrations(system.getComponents().size()) {
				for(auto& component : system.getComponents()) {
					init_concentrations[component->getIndex()] = component->concentration();
				}
			}

		X F::concentrations(const X& extent) const {
			// Init concentrations of all components
			std::vector<double> concentrations = init_concentrations;
			for(auto& reaction : system.getReactions()) {
				const std::vector<double>& coefficients
					= system.getReactionMatrix()[reaction->getIndex()];
				for(
						std::size_t component_index = 0;
						component_index < coefficients.size();
						component_index++) {

					concentrations[component_index] =
						system.getComponent(component_index)
								.concentration(
									concentrations[component_index],
									coefficients[component_index]*extent[reaction->getIndex()]
									);
				}
			}
			// Activities of all components resulting from the extent of all
			// reactions
			return concentrations;
		}
		X F::f(const X& extent) const {
			X concentrations = this->concentrations(extent);

			X f(extent.size());
			for(auto& reaction : system.getReactions()) {
				const std::vector<double>& coefficients
					= system.getReactionMatrix()[reaction->getIndex()];
				for(
						std::size_t component_index = 0;
						component_index < coefficients.size();
						component_index++) {
					auto& component = system.getComponent(component_index);
					f[reaction->getIndex()] +=
						// Stoichiometric coefficient
						coefficients[component_index] * log(
								// Current activity of the component according
								// to the provided reaction extent
								component.activity(concentrations[component_index])
								);
				}
				// The previous sum should be equal to log(K)
				f[reaction->getIndex()] += reaction->getLogK();
			}
			return f;
		}
		M F::df(const X& extent) const {
			X current_concentrations = this->concentrations(extent);
			M jacobian(extent.size()); // reactions count

			for(auto& reaction : system.getReactions()) {
				jacobian[reaction->getIndex()].resize(extent.size());
				const std::vector<double>& coefficients_i
					= system.getReactionMatrix()[reaction->getIndex()];
				for(auto& dreaction : system.getReactions()) {
					const std::vector<double>& coefficients_j
						= system.getReactionMatrix()[dreaction->getIndex()];
					double dfi_dxj = 0;
					for(const auto& component : system.getComponents()) {
						dfi_dxj += coefficients_i[component->getIndex()] *
							component->Dactivity(
									current_concentrations,
									coefficients_j[component->getIndex()]) /
							component->activity(
									current_concentrations[component->getIndex()]
									) / ln10;
					}
					jacobian[reaction->getIndex()][dreaction->getIndex()] = dfi_dxj;
				}
			}

			return jacobian;
		}

		X solve(const ChemicalSystem& system) {
			F f(system);
			// Initially all to 0
			X init_extent(system.getReactions().size());
			for(const auto& reaction : system.getReactions()) {
				double reactives_activity = 1.0;
				double products_activity = 1.0;
				std::cout << "Reaction " << reaction->getName() << std::endl;
				for(const auto& reactive : reaction->getReagents()) {
					std::cout << "  " << reactive.first << "(" << reactive.second << "): " << 
							system.getComponent(reactive.first).activity() << std::endl;
					if(reactive.second < 0.0)
						products_activity *=
							system.getComponent(reactive.first).activity();
					else if(reactive.second > 0.0)
						reactives_activity *=
							system.getComponent(reactive.first).activity();
				}
				if(reactives_activity == 0.0 && products_activity == 0.0) {
					throw std::logic_error("Invalid reaction with missing reactives and products");
				} else {
					if(reactives_activity == 0.0) {
						init_extent[reaction->getIndex()] = 1e-16;
					} else if (products_activity == 0.0) {
						init_extent[reaction->getIndex()] = -1e-16;
					} else {
						init_extent[reaction->getIndex()] = 0.0;
					}
				}
			}
			return Newton<X, M>(
					init_extent,
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(system.getMaxIteration());
		}
	}
	}


	/*
	 *ChemicalSystem::ChemicalSystem(
	 *        double _K1, double _K2, double _K3, double initPH,
	 *        double H, double P, double C,
	 *        double N,
	 *        double S, double SH, double SP, double SC
	 *        ) :
	 *    _K1(_K1), _K2(_K2), _K3(_K3), _initPH(initPH),
	 *    H(H), P(P), C(C), N(N), S(S), SH(SH), SP(SP), SC(SC) {
	 *        std::cout << "[Init] Mineral weight: " << mineral_weight/kg << "kg" << std::endl;
	 *        std::cout << "[Init] H: " << H/mol << "mol = " << (H/mol)/(V/l) << "mol/L" << std::endl;
	 *        std::cout << "[Init] P: " << P/mol << "mol = " << (P/mol)/(V/l) << "mol/L = " <<
	 *            ((P/entities)*(30.973*u)/mg) / (mineral_weight/kg) << "mg/kg" << std::endl;
	 *        std::cout << "[Init] C: " << C/mol << "mol = " << (C/mol)/(V/l) << "mol/L = " <<
	 *            ((C/entities)*(12.010*u)/mg) / (mineral_weight/kg) << "mg/kg" << std::endl;
	 *        std::cout << "[Init] pH: " << this->initPH() << std::endl;
	 *        std::cout << std::endl;
	 *        std::cout << "[Init] N: " << N << " sites (mol)." << std::endl;
	 *        std::cout << "[Init] S : " << S/mol << " sites (mol)." << std::endl;
	 *        std::cout << "[Init] SH: " << SH/mol << " sites (mol)." << std::endl;
	 *        std::cout << "[Init] SP: " << SP/mol << " sites (mol)." << std::endl;
	 *        std::cout << "[Init] SC: " << SC/mol << " sites (mol)." << std::endl;
	 *    }
	 */

/*
 *    ChemicalSystem ChemicalSystem::equilibrium(
 *            double K1, double K2, double K3, double pH,
 *            double solution_P, double solution_C,
 *            double mineral_N
 *            ) {
 *        double _K1 = K1/V;
 *        double _K2 = K2/std::pow(V, 2);
 *        double _K3 = K3/V;
 *        double H = V*std::pow(10, -pH)*mol/l;
 *        double P = ((mineral_weight * solution_P)/(30.973*u)) * entities;
 *        double C = ((mineral_weight * solution_C)/(12.010*u)) * entities;
 *        double N = mineral_N;
 *
 *        // Assumes that the system starts at equilibrium
 *        double S = N/(1+_K1*H+_K2*H*P+_K3*C); // S+SH+SP+SC=N
 *        double SH = _K1 * S * H;
 *        double SP = _K2 * S * P * H;
 *        double SC = _K3 * S * C;
 *
 *        return ChemicalSystem(
 *                _K1, _K2, _K3, pH,
 *                H, P, C, N,
 *                S, SH, SP, SC);
 *    }
 */

	/*
	 *ChemicalSystem ChemicalSystem::defaultEquilibrium() {
	 *    return equilibrium(
	 *            std::pow(10, 8.5), // K1
	 *            std::pow(10, 26.3), // K2
	 *            std::pow(10, 0.6), // K3
	 *            7.5, // pH
	 *            1.43e-6 * gram / gram, // 1.43e-6 grams of P in solution by gram of soil
	 *            729.0e-6 * gram / gram, // 729.0e-6 grams of C in solution by gram of soil
	 *            1.45e19 * entities / gram * mineral_weight // 1.45e19 sites by gram of soil
	 *            );
	 *}
	 */

	/*
	 *ChemicalSystem ChemicalSystem::Devau2011Control() {
	 *    return equilibrium(
	 *            // Kaolinite
	 *            std::pow(10, 4.36),
	 *            std::pow(10, 23),
	 *            std::pow(10, 1),
	 *            6.5,
	 *            (279-69)*mg/kg,
	 *            9.80*gram/kg,
	 *            (6.15*entities/nm2)*(105*m2/gram)*(20.12*gram/l)*V
	 *            );
	 *}
	 */

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

