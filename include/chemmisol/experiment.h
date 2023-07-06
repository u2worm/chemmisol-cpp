#include "chemical.h"
#include <fstream>
#include <functional>
#include <unordered_map>

namespace chemmisol {
	class Experiment {
		// TODO: experiment parameters

		const ChemicalSystem& base_system;
		std::ofstream csv_file;
		std::function<double(unsigned long)> pH;
		unsigned long t = 0;

		private:
		void writeToCsv(
				const std::unordered_map<std::string, std::vector<double>>& results);

		public:
		Experiment(
				std::string output_file,
				const ChemicalSystem& system,
				const std::function<double(unsigned long)>& pH
				);
		std::unordered_map<std::string, std::vector<double>> run(unsigned long time);
	};
}
