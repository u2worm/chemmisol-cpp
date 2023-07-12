#ifndef CHEMMISOL_LOGGING_H
#define CHEMMISOL_LOGGING_H

#include "easylogging++.h"

#define CHEM_CORE_LOGID "chemmisol-core"
#define CHEM_LOG(LEVEL) CLOG(LEVEL, CHEM_CORE_LOGID)
#define CHEM_LOGV(VERBOSE_LEVEL) CVLOG(VERBOSE_LEVEL, CHEM_CORE_LOGID)

namespace chemmisol {
	extern el::Logger* chemmisol_core_logger;
	extern el::Configurations chemmisol_core_logger_configuration;

	el::Logger* init_logger(const char* logger_id);
	void log_file(el::Logger* logger, const std::string& _log_file, bool terminal_logging);
	void log_file(const std::string& file, bool terminal_logging = false);
}
#endif
