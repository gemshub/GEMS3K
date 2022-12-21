#include <set>
#include <spdlog/sinks/stdout_color_sinks.h>
#include "v_detail.h"

// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> gems_logger = spdlog::stdout_color_mt("gems3k");

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

