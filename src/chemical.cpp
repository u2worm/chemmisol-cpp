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

	bool operator==(const ReactionComponent& c1, const ReactionComponent& c2) {
		return c1.name == c2.name
			&& c1.phase == c2.phase
			&& c1.coefficient == c2.coefficient;
	}

	ChemicalSystem::ChemicalSystem(
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					const std::string& surface_component
					) :
		solid_concentration(solid_concentration),
		specific_surface_area(specific_surface_area),
		site_concentration(site_concentration) {
			addComponent(surface_component, MINERAL, 1);
		}

	ChemicalSystem::ChemicalSystem(const ChemicalSystem& other) :
				// Indexes initialized to 0 since components and reactions are
				// copied using addComponent() and addReaction()
				component_index(0),
				reaction_index(0),
				max_iteration(other.max_iteration),
				solid_concentration(other.solid_concentration),
				specific_surface_area(other.specific_surface_area),
				site_concentration(other.site_concentration) {
					for(auto& component : other.getComponents()) {
						if(component->isFixed()) {
							this->fixComponent(
									component->getName(),
									component->getPhase(),
									component->concentration()
									);
						} else {
							this->addComponent(
									component->getName(),
									component->getPhase(),
									component->concentration()
									);
						}
					}
					for(auto& reaction : other.getReactions()) {
						this->addReaction(
								reaction->getName(),
								reaction->getLogK(),
								reaction->getReagents()
								);
					}
					this->reaction_matrix = other.reaction_matrix;
			}

	void ChemicalSystem::addReaction(Reaction* reaction) {
		reactions.emplace_back(reaction);
		reactions_by_name.emplace(reaction->getName(), reaction);
	}

	void ChemicalSystem::addReaction(
			std::string name,
			double K,
			std::vector<ReactionComponent> reactives
			) {
		addReaction(new Reaction(
				name, reaction_index++, K, reactives
				));
	}

	void ChemicalSystem::addComponent(Component* component) {
		components.emplace_back(component);
		components_by_name.emplace(component->getName(), component);
	}

	void ChemicalSystem::fixComponent(
			const std::string& name,
			double concentration
			) {
		fixComponent(name, AQUEOUS, concentration);
	}

	void ChemicalSystem::fixComponent(
			const std::string& name,
			Phase phase,
			double concentration
			) {

		auto existing_component = components_by_name.find(name);
		if (existing_component != components_by_name.end()) {
			std::size_t index = existing_component->second->getIndex();
			FixedComponent* fixed_component;
			switch(phase) {
				case AQUEOUS:
					fixed_component = new FixedAqueousComponent(
							name, index,
							concentration
							);
					break;
				case MINERAL:
					fixed_component = new FixedMineralComponent(
							name, index,
							concentration
							);
					break;
				case SOLVENT:
					fixed_component = new Solvent(name, index);
			}
			components[index].reset(fixed_component);
			existing_component->second = fixed_component;
		} else {
			switch(phase) {
				case AQUEOUS:
					addComponent(new FixedAqueousComponent(
							name, component_index++,
							concentration
							));
					break;
				case MINERAL:
					addComponent(new FixedMineralComponent(
							name, component_index++,
							concentration
							));
					break;
				case SOLVENT:
					addComponent(new Solvent(name, component_index++));
			}
		}
	}

	void ChemicalSystem::addComponent(const std::string& name, double concentration) {
		addComponent(name, AQUEOUS, concentration);
	}

	void ChemicalSystem::addComponent(
			const std::string& name, Phase phase, double concentration) {
		switch(phase) {
			case AQUEOUS:
				addComponent(
						new AqueousComponent(name, component_index++, concentration)
						);
				break;
			case MINERAL:
				addComponent(
						new MineralComponent(
							name, component_index++,
							solid_concentration, specific_surface_area, site_concentration,
							concentration)
						);
				break;
			case SOLVENT:
				addComponent(new Solvent(name, component_index++));
		}
	}

	void ChemicalSystem::initPH(double pH) {
		addComponent("H+", std::pow(10, -pH)*mol/l);
		addComponent("OH-", std::pow(10, -pH)*mol/l);
	}
	void ChemicalSystem::fixPH(double pH) {
		fixComponent(
					"H+", AQUEOUS,
					std::pow(10, -pH)*mol/l
				);
		addComponent(
					"OH-", AQUEOUS,
					std::pow(10, getReaction("OH-").getLogK())
						/ std::pow(10, -pH)*mol/l
				);
		//initPH(pH);
		// A fake reaction that produces (or consumes) as much H+ as needed to
		// fix the pH to the provided value
		//addReaction("pH", -pH, {{"H+", -1}});
	}

	double ChemicalSystem::getPH() const {
		return -log(getComponent("H+").concentration());
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
				if(components_by_name.find(component.name) == components_by_name.end()) {
					addComponent(component.name, component.phase, 0);
				}
			}
		}

		reaction_matrix.resize(reactions.size());
		for(const auto& reaction : reactions) {
			reaction_matrix[reaction->getIndex()].resize(components.size());
			for(const auto& reactive : reaction->getReagents()) {
				reaction_matrix[reaction->getIndex()][getComponent(reactive.name).getIndex()]
					= reactive.coefficient;
			}
		}
	}

	void ChemicalSystem::solveEquilibrium() {
		initReactionMatrix();
		using namespace solver;
		solver::X dx = solve(*this);
		for(auto& reaction : reactions) {
			for(auto& component : reaction->getReagents()) {
				components[getComponent(component.name).getIndex()]->incrementConcentration(
						component.coefficient * dx[reaction->getIndex()]
						);
			}
		}
	}

	double ChemicalSystem::reactionQuotient(const std::string& name) const {
		double quotient = 1;
		auto& reaction = getReaction(name);
		for(auto& reactive : reaction.getReagents()) {
			quotient *= std::pow(getComponent(reactive.name).activity(), reactive.coefficient);
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
					std::cout << "  " << reactive.name << "(" << reactive.coefficient << "): " << 
							system.getComponent(reactive.name).activity() << std::endl;
					if(reactive.coefficient < 0.0)
						products_activity *=
							system.getComponent(reactive.name).activity();
					else if(reactive.coefficient > 0.0)
						reactives_activity *=
							system.getComponent(reactive.name).activity();
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

