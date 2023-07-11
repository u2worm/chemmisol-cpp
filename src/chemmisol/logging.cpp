#include "chemmisol/logging.h"

INITIALIZE_EASYLOGGINGPP

namespace chemmisol {
	static el::Logger* init_logger();
	void log_file(el::Logger*, const std::string&, bool);

	el::Logger* chemmisol_core_logger = init_logger();

	static el::Logger* init_logger() {
		el::Logger* logger = el::Loggers::getLogger(CHEM_CORE_LOGID);

		el::Configurations chemmisol_core_logger_configuration;
		// Global format
		chemmisol_core_logger_configuration.setGlobally(
				el::ConfigurationType::Format, "[%logger] %levshort %msg");

		//defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "false");

		el::Loggers::reconfigureLogger(logger, chemmisol_core_logger_configuration);

		//el::Loggers::setVerboseLevel(2);
   		return logger;
	}

	void log_file(el::Logger* logger, const std::string& _log_file, bool terminal_logging) {
		el::Configurations chemmisol_core_logger_configuration;
		chemmisol_core_logger_configuration.setGlobally(
				el::ConfigurationType::ToFile, "true");
		chemmisol_core_logger_configuration.setGlobally(
				el::ConfigurationType::Filename, _log_file);

		if(!terminal_logging)
			chemmisol_core_logger_configuration.setGlobally(
					el::ConfigurationType::ToStandardOutput, "false");

		el::Loggers::reconfigureLogger(
				logger, chemmisol_core_logger_configuration);
	}


	void log_file(const std::string& _log_file, bool disable_terminal_logging) {
		log_file(chemmisol_core_logger, _log_file, disable_terminal_logging);
	}
}
