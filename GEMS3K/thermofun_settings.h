// The settings of a ThermoFun engine that works for GEMS.
#ifndef THERMOFUN_SETTINGS_H
#define THERMOFUN_SETTINGS_H

#ifndef NO_THERMOFUN

#include "ThermoFun/ThermoFun.h"

namespace thermofun_settings
{

// Older versions of ThermoFun do not have this option: their G has no fugacity term.
template <class Preferences>
auto cancelFugacity(Preferences& p, int) -> decltype(p.cancel_fugacity_coefficient_in_gas_props = true, void())
{
    p.cancel_fugacity_coefficient_in_gas_props = true;
}
template <class Preferences>
void cancelFugacity(Preferences&, long) {}

// GEMS adds the pressure and the fugacity terms of a gas itself, so the engine gives G(T, Pr) of a gas. From ThermoFun 1.0 both terms are on by default.
inline void useStandardStateOfGases(ThermoFun::ThermoEngine& engine)
{
    // The reference is changed in place: a libThermoFun older than its installed header has no setPreferences().
    auto& preferences = engine.preferences();
    preferences.apply_pressure_correction_to_gas_props = false;
    cancelFugacity(preferences, 0);
}

} // namespace thermofun_settings

#endif
#endif
