#include "pySolmod4rkt.hpp"

PYBIND11_MODULE(PySolmod4rkt, m)
{
    // Read config file
    gemsSettings().gems3k_update_loggers(false, "solmod.log", 3);
    exportSolMod(m);
}
