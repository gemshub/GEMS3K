#include <set>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include "v_detail.h"

// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> gems_logger = spdlog::stdout_color_mt("gems3k");

std::set<std::string> gems3k_loggers = {
    "gems3k", "ipm", "tnode", "kinmet", "solmod",
#ifdef USE_THERMOFUN
    "chemicalfun", "thermofun"
#endif
};
std::string gems3k_logger_pattern("[%n] [%^%l%$] %v");

void gems3k_update_loggers( bool use_cout, const std::string& logfile_name, size_t log_level)
{
    std::set<std::shared_ptr<spdlog::logger>> gems3k_loggers_set;
    for(const auto& lname: gems3k_loggers) {
        auto res = spdlog::get(lname);
        if (res) {
            gems3k_loggers_set.insert(res);
        }
    }

    spdlog::level::level_enum log_lev = spdlog::level::info;
    if( log_level<7 ) {
        log_lev = static_cast<spdlog::level::level_enum>(log_level);
    }

    for(const auto& logger: gems3k_loggers_set) {
        logger->sinks().clear();
        if(use_cout) {
            auto console_output = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
            console_output->set_pattern(gems3k_logger_pattern);
            logger->sinks().push_back(console_output);
        }
        if(!logfile_name.empty()) {
            auto file_output = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile_name, 1048576, 3);
            file_output->set_pattern(gems3k_logger_pattern);
            logger->sinks().push_back(file_output);
        }
        logger->set_level(log_lev);
    }
}

void gems3k_clear_loggers( const std::string& logfile_name)
{
    std::set<std::shared_ptr<spdlog::logger>> gems3k_loggers_set;
    for(const auto& lname: gems3k_loggers) {
        auto res = spdlog::get(lname);
        if (res) {
            gems3k_loggers_set.insert(res);
        }
    }

    for(const auto& logger: gems3k_loggers_set) {
        logger->sinks().clear();
        if(!logfile_name.empty()) {
            auto file_output = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile_name, 1048576, 3);
            logger->sinks().push_back(file_output);
        }
    }
}

TError::~TError()
{}

[[ noreturn ]] void Error (const std::string& title, const std::string& message)
{
    gems_logger->error("{}: {}", title, message);
    throw TError(title, message);
}

void ErrorIf (bool error, const std::string& title, const std::string& message)
{
    if(error) {
        gems_logger->error("{}: {}", title, message);
        throw TError(title, message);
    }
}

template <> double InfMinus()
{
  return DOUBLE_INFMINUS;
}
template <> double InfPlus()
{
  return DOUBLE_INFPLUS;
}
template <> double Nan()
{
  return DOUBLE_NAN;
}

template <> float InfMinus()
{
  return FLOAT_INFMINUS;
}
template <> float InfPlus()
{
  return FLOAT_INFPLUS;
}
template <> float Nan()
{
  return FLOAT_NAN;
}

