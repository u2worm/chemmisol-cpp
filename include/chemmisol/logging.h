#ifndef CHEMMISOL_LOGGING_H
#define CHEMMISOL_LOGGING_H

#include "easylogging++.h"

#define CHEM_CORE_LOGID "chemmisol-core"
#define CHEM_LOG(LEVEL) CLOG(LEVEL, CHEM_CORE_LOGID)
#define CHEM_LOGV(VERBOSE_LEVEL) CVLOG(VERBOSE_LEVEL, CHEM_CORE_LOGID)

namespace chemmisol {
	extern el::Logger* chemmisol_core_logger;
	extern el::Configurations chemmisol_core_logger_configuration;

	void log_file(const std::string& file, bool terminal_logging = false);
}
#endif
