#include "experiment.h"
#include <unordered_map>

namespace mineral {
	Experiment::Experiment(
			std::string output_file,
			const ChemicalSystem& chemical_system,
			const std::function<double(unsigned long)>& pH
			)
		: base_system(chemical_system), csv_file(output_file), pH(pH) {
		}

	void Experiment::writeToCsv(
			const std::unordered_map<std::string, std::vector<double>>& results) {
		csv_file << "pH";
		for(auto& item : results)
			csv_file << "," << item.first;
		csv_file << std::endl;

		for(std::size_t i = 0; i < results.begin()->second.size(); i++) {
			csv_file << -std::log10(results.at("H+")[i]);
			for(auto& item : results) {
				csv_file << "," << item.second[i];
			}
			csv_file << std::endl;
		}
	}

	std::unordered_map<std::string, std::vector<double>> Experiment::run(unsigned long time) {
		std::unordered_map<std::string, std::vector<double>> results(time);

		t = 0;
		while(t < time) {
			ChemicalSystem system_copy(base_system);
			std::cout << "Exp pH: " << pH(t) << std::endl;
			system_copy.fixPH(pH(t));

			system_copy.solveEquilibrium();

			for(auto& component : system_copy.getComponents())
				results[component->getName()].push_back(
						component->concentration()
						);
			++t;
		}
		writeToCsv(results);
		return results;
	}
}
