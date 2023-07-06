#include "chemmisol/chemical.h"
#include "chemmisol/newton.h"
#include "chemmisol/regula_falsi.h"

#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>

namespace chemmisol {
	/*
	 *const double ChemicalSystem::bulk_density = 1.17*gram/cm3;
	 *const double ChemicalSystem::V = std::pow(1*cm, 3);
	 *const double ChemicalSystem::mineral_weight  = V*bulk_density;
	 */

	const double Component::V = 1*l;

	double MineralComponent::sites_count(
			double solid_concentration,
			double specific_surface_area,
			double site_concentration) {
		return V*solid_concentration*specific_surface_area*site_concentration;
	}

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
							solid_concentration, specific_surface_area,
							site_concentration,
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
							solid_concentration, specific_surface_area,
							site_concentration,
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
		//addComponent(
					//"OH-", AQUEOUS,
					//std::pow(10, getReaction("OH-").getLogK())
						/// std::pow(10, -pH)*mol/l
				//);
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

	double GuessF::operator()(const double &extent) const {
		double f = 0;
		for(auto& reagent : reaction.getReagents()) {
			//std::cout << "    " << reagent.name << " cc: " << current_concentrations[system.getComponent(reagent.name).getIndex()] << " + " << reagent.coefficient * extent
				//<< " = " << system.getComponent(reagent.name).activity(
					//current_concentrations,
					//reagent.coefficient * extent
					//)
				//<< std::endl;
			f += reagent.coefficient * log(system.getComponent(reagent.name).activity(
					current_concentrations,
					reagent.coefficient * extent
					));
		}
		f += reaction.getLogK();
		return f;
	}

	const std::unordered_map<std::string, double>& ChemicalSystem::guessInitialExtents() {
		for(const auto& reaction : getReactions()) {
			if(initial_guess_extents.find(reaction->getName())
					== initial_guess_extents.end()) {
				// No manually specified guess extent: automatically compute one

				double reactants_product = 1.0;
				for(const auto& reagent : reaction->getReagents()) {
					if(reagent.coefficient > 0.0) {
						reactants_product *=
							std::pow(getComponent(reagent.name).activity(), reagent.coefficient);
					}
				}
				if(reactants_product == 0.0) {
					// TODO: support for this situation?
					throw std::logic_error(
							"Some reactants are missing: the problem is ill formed."
							);
				}
			}
		}
		// Concentrations used in the incremental guess process
		std::vector<double> guessed_concentrations(getComponents().size());
		for(auto& component : getComponents())
			guessed_concentrations[component->getIndex()]
				= component->concentration();

		// The overall idea of the guess process is that we have an initial
		// quantity of reactives, specified by the user, and some null products.
		// We need to:
		// 1. initialize the products quantity to a non null value.
		// 2. initialize it to a value as close as possible to the final
		// equilibrium so the Newton method can converge.
		//
		// The logK() values are the first hint to find activities. Considering
		// that all free products concentrations are null, reactions
		// with a big positive logK() value will try to consume as much
		// reactive as possible to maximise the products quantity and establish
		// the equilibrium, while reactions with a big negative logK() value
		// will only need to consumme a negligible quantity of reactives to
		// establish the equilibrium.
		//
		// In order to estimate the final equilibrium, we start by ordering
		// reactions by their logK() value, so that reactions that need a lot of
		// reactives will be handled first. Then the algorithm consists in
		// solving each reaction as if they where alone in the system. At each
		// reaction, the estimated quantity of reactives consumed is updated in
		// guessed_concentrations, so that the first reaction with big logK()
		// value have a lot of reactives at the beginning, and the last one can
		// still get close to the equilibrium from what's left.

		std::vector<const Reaction*> sorted_reactions;
		{
			// Sort reactions by logK() values
			for(auto& reaction : getReactions())
				sorted_reactions.push_back(reaction.get());

			std::sort(sorted_reactions.begin(), sorted_reactions.end(),
					[] (const Reaction*& r1, const Reaction*& r2) {
					// > instead of < so that the list is sorted in descending
					// order
					return r1->getLogK() > r2->getLogK();
					});
		}

		for(auto& reaction : sorted_reactions) {
			std::cout <<
				"Try to guess extent for reaction " << reaction->getName() <<
				" (log K = " << reaction->getLogK() << ")" << std::endl;

			// Limiting reactive finding algorithm
			//
			// The limiting reactive is used to define the maximum possible
			// extent of the reaction so that all reagent concentrations stay
			// positive.
			const Component* limiting_reactive;
			double limiting_reactive_coefficient = 0;
			double smallest_limiting_factor = std::numeric_limits<double>::infinity();

			for(auto& reagent : reaction->getReagents()) {
				if(reagent.coefficient > 0.0) {
					// Consider only reactives
					const Component& component = getComponent(reagent.name);
					if(!component.isFixed()) {
						// A fixed reactive cannot be limiting
						double limiting_factor = component.quantity(
								guessed_concentrations[component.getIndex()]
								) / reagent.coefficient;
						//std::cout << "  v " << reagent.coefficient << " c " << component.quantity(guessed_concentrations[component.getIndex()]) << " " << component.getName() << " ~ " << limiting_factor << std::endl;
						if(limiting_factor < smallest_limiting_factor) {
							limiting_reactive = &component;
							limiting_reactive_coefficient = reagent.coefficient;
							smallest_limiting_factor = limiting_factor;
						}
					}
				} 
			}

			std::cout << "  Limiting reactive: " << limiting_reactive_coefficient
				<< " " << limiting_reactive->getName() << " " << 
				limiting_reactive_coefficient * smallest_limiting_factor << " mol/l"
				<< std::endl;

			// See the GuessF definition for more details.
			// Return log(Q(x)) - logK where x is the extent of the reaction and
			// Q is the reaction quotient with products at the numerator
			GuessF f = GuessF(*this, *reaction, guessed_concentrations);

			double max_N;
			if(smallest_limiting_factor < std::numeric_limits<double>::infinity()) {
				max_N = smallest_limiting_factor;
			} else {
				// In this case, there is no limiting factor. This for example
				// the case for the reaction OH-: H2O <-> H+ + OH- where H2O is
				// a solvent and H+ is fixed with the pH. Using a value of max_N
				// = 1mol/l still allows to find the corresponding OH-
				// concentration using the RegulaFalsi solver.
				max_N = 1*mol/l;
			}

			//std::cout << 
					//-(max_N) * (1 - std::numeric_limits<double>::epsilon())
					//<< std::endl;

			// Uses the RegulaFalsi method to find the root of f(x)=log(Q(x)) - logK
			// Solutions are searched in ]-N, 0[ where N is the limiting factor.
			// Due to the convention that products are specified with negative
			// coefficients and that we find a solution that necessarily form
			// products, extents considered are negative.
			double guess_extent = RegulaFalsi<double>(
					// Value just below 0
					- std::numeric_limits<double>::min(),
					// Value just above -N
					//-(max_N) * (1 - std::numeric_limits<double>::epsilon()),
					-(max_N) * (1 - std::numeric_limits<double>::epsilon()),
					// Function to optimize
					f
					).solve_iter(10000);
			// About std::numeric_limits<double>::epsilon():
			// epsilon is the next float that can be represented after 1.0, so
			// that it is the smallest quantity such that:
			// 1.0 - epsilon < 1.0 < 1.0 + epsilon
			// Multiplying this inequality by N gives the "just above -N" value.

			std::cout << "  Guessed extent: " << guess_extent << std::endl;

			// Set up the guessed value
			initial_guess_extents[reaction->getName()] = guess_extent;
			// Updates concentrations from the guessed extent, so that the
			// system remains consistent for next reactions.
			for(auto& reagent : reaction->getReagents()) {
				auto& component = getComponent(reagent.name);
				guessed_concentrations[component.getIndex()]
					= component.concentration(
							guessed_concentrations[component.getIndex()],
							reagent.coefficient * guess_extent
							);
			}
		}
		return initial_guess_extents;
	}

	void ChemicalSystem::solveEquilibrium() {
		initReactionMatrix();
		guessInitialExtents();

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
				std::cout << "Calc C from reaction " << reaction->getName() << std::endl;
				for(
						std::size_t component_index = 0;
						component_index < coefficients.size();
						component_index++) {

					std::cout << "  " <<
						system.getComponent(component_index).getName() << ": " <<
						concentrations[component_index] << " + " <<
						coefficients[component_index]*extent[reaction->getIndex()] << " = " <<
						system.getComponent(component_index)
								.concentration(
									concentrations[component_index],
									coefficients[component_index]*extent[reaction->getIndex()]
									)
						<< std::endl;
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
			std::cout << "F concentrations: " << std::endl;
			for(auto& component : system.getComponents())
				std::cout << "  " << component->getName() << " " << concentrations[component->getIndex()] << std::endl;
			std::cout << std::endl;

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
				/*
				 *double reactives_activity = 1.0;
				 *double products_activity = 1.0;
				 *std::cout << "Reaction " << reaction->getName() << std::endl;
				 *for(const auto& reactive : reaction->getReagents()) {
				 *    std::cout << "  " << reactive.name << "(" << reactive.coefficient << "): " << 
				 *            system.getComponent(reactive.name).activity() << std::endl;
				 *    if(reactive.coefficient < 0.0)
				 *        products_activity *=
				 *            system.getComponent(reactive.name).activity();
				 *    else if(reactive.coefficient > 0.0)
				 *        reactives_activity *=
				 *            system.getComponent(reactive.name).activity();
				 *}
				 *if(reactives_activity == 0.0 && products_activity == 0.0) {
				 *    throw std::logic_error("Invalid reaction with missing reactives and products");
				 *} else {
				 *    if(reactives_activity == 0.0) {
				 *        init_extent[reaction->getIndex()] = 1e-8;
				 *    } else if (products_activity == 0.0) {
				 *        init_extent[reaction->getIndex()] = -1e-8;
				 *    } else {
				 *        init_extent[reaction->getIndex()] = 0.0;
				 *    }
				 *}
				 */
				init_extent[reaction->getIndex()]
					= system.getInitialGuessExtent(reaction->getName());
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

