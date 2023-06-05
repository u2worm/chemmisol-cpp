#include "chemical.h"
#include <fstream>
#include <functional>

namespace mineral {
	class Experiment {
		// TODO: experiment parameters

		ChemicalSystem system;
		std::ofstream csv_file;
		std::function<double(unsigned long)> pH;
		std::function<double(unsigned long, double)> dP;
		unsigned long t = 0;

		private:
		void writeToCsv();

		public:
		Experiment(
				std::string output_file,
				const ChemicalSystem& system,
				const std::function<double(unsigned long)>& pH,
				const std::function<double(unsigned long, double)>& dP
				);
		void run(unsigned long time);
	};
}
