#include "chemmisol/logging.h"

INITIALIZE_EASYLOGGINGPP

namespace chemmisol {
	void log_file(el::Logger*, const std::string&, bool);

	el::Logger* core_logger = init_logger(CHEM_CORE_LOGID);

	el::Logger* init_logger(const char* logger_id) {
		el::Logger* logger = el::Loggers::getLogger(logger_id);

		el::Configurations core_logger_configuration;
		// Global format
		core_logger_configuration.setGlobally(
				el::ConfigurationType::Format, "[%logger] %levshort %msg");

		//defaultConf.set(el::Level::Trace, el::ConfigurationType::Enabled, "false");

		el::Loggers::reconfigureLogger(logger, core_logger_configuration);

		//el::Loggers::setVerboseLevel(2);
   		return logger;
	}

	void log_to_file(el::Logger* logger, const std::string& _log_file, bool terminal_logging) {
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
}
