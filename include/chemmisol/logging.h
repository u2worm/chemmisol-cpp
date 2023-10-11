#ifndef CHEMMISOL_LOGGING_H
#define CHEMMISOL_LOGGING_H

#include "easylogging++.h"

/**
 * @file chemmisol/logging.h
 *
 * Logging features resulting from the integration of
 * [easylogging++](https://github.com/abumq/easyloggingpp).
 */

/**
 * Logging ID used by the chemmisol::core_logger.
 */
#define CHEM_CORE_LOGID "chemmisol-core"
/**
 * [Easylogging++ LOG](https://github.com/abumq/easyloggingpp#logging) to the
 * chemmisol::core_logger.
 */
#define CHEM_LOG(LEVEL) \
	CLOG(LEVEL, CHEM_CORE_LOGID)

/**
 * [Easylogging++ VLOG](https://github.com/abumq/easyloggingpp#verbose-logging)
 * to the chemmisol::core_logger.
 */
#define CHEM_LOGV(VERBOSE_LEVEL) CVLOG(VERBOSE_LEVEL, CHEM_CORE_LOGID)

namespace chemmisol {
	/**
	 * Pointer to the Easylogging++ Logger used internally by the chemmisol
	 * library.
	 */
	extern el::Logger* core_logger;

	/**
	 * Initializes and configures a chemmisol logger with the specified ID.
	 *
	 * The method is notably used to initialize the chemmisol::core_logger with
	 * the #CHEM_CORE_LOGID, but can be used to define other loggers with the
	 * same configuration ad the chemmisol::core_logger.
	 *
	 * @param logger_id ID of the logger, used in Easylogging++ CLOG() and
	 * CVLOG() macros.
	 * @return Initialized and configured logger.
	 */
	el::Logger* init_logger(const char* logger_id);

	/**
	 * Configures the specified logger to output to the specified log file.
	 * Logging to std::cout can also be disabled thanks to the terminal_logging
	 * argument.
	 *
	 * Any Easylogging++ logger can be configured this way, including the
	 * chemmisol::core_logger or any logger initialized with init_logger().
	 *
	 * @param logger Logger to output to a file.
	 * @param log_file Full path to the log file.
	 * @param terminal_logging Set to false to disable std::cout logging.
	 */
	void log_to_file(
			el::Logger* logger, const std::string& log_file, bool terminal_logging
			);
}
#endif
