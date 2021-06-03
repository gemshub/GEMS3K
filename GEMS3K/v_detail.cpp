#include "v_detail.h"

TError::~TError()
{}

[[ noreturn ]] void Error (const std::string& title, const std::string& message)
{
    throw TError(title, message);
}

void ErrorIf (bool error, const std::string& title, const std::string& message)
{
    if(error)
        throw TError(title, message);
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

