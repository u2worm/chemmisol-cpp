#include "chemmisol/chemical/reaction.h"
#include <sstream>
#include "chemmisol/chemical/system.h"

namespace chemmisol {

	EmptyReagents::EmptyReagents(
			const ChemicalSystem* chemical_system,
			const Reaction* invalid_reaction) :
		InvalidReaction(chemical_system, invalid_reaction) {
			std::ostringstream message_stream;
			message_stream << "Empty reagents list in the reaction " <<
				invalid_reaction->getName() << "."
				<< std::endl;
			message_stream << "Components currently declared in the chemical system :"
				<< std::endl;
			for(const auto& component : chemical_system->getComponents()) {
				message_stream << "  "
					<< component->getSpecies()->getName() << " "
					<< "(" << component->getSpecies()->getPhase() << ")"
					<< std::endl;
			}
			message = message_stream.str();
		}
	MissingProducedSpeciesInReaction::MissingProducedSpeciesInReaction(
			const ChemicalSystem* chemical_system,
			const Reaction* invalid_reaction) :
		InvalidReaction(chemical_system, invalid_reaction) {
			std::ostringstream message_stream;
			message_stream << "No produced species could be found for the "
				"reaction " << invalid_reaction->getName() << ", and components "
				"only reactions are not allowed. Check your reaction definition "
				"and ensure that no extra species are declared as components."
				<< std::endl;
			message_stream << "Reagents currently declared in the reaction "
				<< invalid_reaction->getName() << " :" << std::endl;
			for(const auto& reagent : invalid_reaction->getReagents())
				message_stream << "  " << reagent.coefficient << " " << reagent.name
					<< " (" << reagent.phase << ")" << std::endl;
			message_stream << "Components currently declared in the chemical system :"
				<< std::endl;
			for(const auto& component : chemical_system->getComponents()) {
				message_stream << "  "
					<< component->getSpecies()->getName() << " "
					<< "(" << component->getSpecies()->getPhase() << ")"
					<< std::endl;
			}
			message = message_stream.str();
		}

	TooManyProducedSpeciesInReaction::TooManyProducedSpeciesInReaction(
			const ChemicalSystem* chemical_system,
			const Reaction* invalid_reaction) :
		InvalidReaction(chemical_system, invalid_reaction) {
			std::ostringstream message_stream;
			message_stream << "Too many produced species detected in the reaction " <<
				invalid_reaction->getName() << ". Reactions must currently be "
				"specified in a canonical form, i.e. each reaction equation must "
				"include exactly one compound chemical species (the \"produced "
				"species\") and chemical components. Either check your "
				"components definition, or refactor your chemical equation "
				"system." << std::endl;
			std::list<Reagent> compound_reagents;
			for(const auto& reagent : invalid_reaction->getReagents()) {
				if(!chemical_system->isComponent(reagent.name))
					compound_reagents.push_back(reagent);
			}
			message_stream << compound_reagents.size() <<
				" compound species were currently found in the reaction "
				" specification : [";
			auto it = compound_reagents.begin();
			while(it != compound_reagents.end()) {
				message_stream << it->coefficient << " " << it->name;
				++it;
				if(it != compound_reagents.end())
					message_stream << ",";
			}
			message_stream << "]" << std::endl;

			message_stream << "Components currently declared in the chemical system :"
				<< std::endl;
			for(const auto& component : chemical_system->getComponents()) {
				message_stream << "  "
					<< component->getSpecies()->getName() << " "
					<< "(" << component->getSpecies()->getPhase() << ")"
					<< std::endl;
			}
			message = message_stream.str();
		}
}
