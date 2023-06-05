#include "units.h"
#include <argparse/argparse.hpp>

int main(int argc, char *argv[])
{
	argparse::ArgumentParser program("convert");
	program.add_argument("value")
		.help("Value to convert")
		.scan<'g', double>();
	program.add_argument("-f", "--from")
		.help("Source unit");
	program.add_argument("-t", "--to")
		.help("Target unit");

	try {
		program.parse_args(argc, argv);
	}
	catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << program;
		return 1;
	}
	double value = program.get<double>("value");
	std::string from = program.get<std::string>("--from");
	std::string to = program.get<std::string>("--to");
	double converted = mineral::convert(value, from, to);

	std::cout << value << from << " = " << converted << to << std::endl;
	
	return 0;
}
