#include <set>
#include "v_detail.h"
#include <spdlog/sinks/stdout_color_sinks.h>

// Thread-safe logger to stdout with colors
std::shared_ptr<spdlog::logger> gems_logger = spdlog::stdout_color_mt("gems3k");

TError::~TError()
{}

// Titles such as "E04IPM: Mass Balance Refinement: " already end with a colon.
static std::string logTitle( const std::string& title )
{
    const auto end = title.find_last_not_of( ": " );
    return end == std::string::npos ? title : title.substr( 0, end + 1 );
}

[[ noreturn ]] void Error (const std::string& title, const std::string& message)
{
    gems_logger->error("{}: {}", logTitle( title ), message);
    throw TError(title, message);
}

void ErrorIf (bool error, const std::string& title, const std::string& message)
{
    if(error) {
        gems_logger->error("{}: {}", logTitle( title ), message);
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

