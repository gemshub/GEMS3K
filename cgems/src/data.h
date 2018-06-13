#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <map>
#include <sstream>

struct DataInfo {
public:
    std::map<std::string, std::string> ClassCodes;
    std::map<std::string, double> AtomicMasses;

    DataInfo() {
        ClassCodes["C"] = "e";
        ClassCodes["O"] = "o";
        ClassCodes["H"] = "h";
        ClassCodes["Ca"] = "e";
        ClassCodes["Mg"] = "e";

        AtomicMasses["H"] = 	0.00100794;
        AtomicMasses["He"] = 	0.004002602;
        AtomicMasses["Li"] = 	0.006941;
        AtomicMasses["Be"] = 	0.009012182;
        AtomicMasses["B"] = 	0.010811;
        AtomicMasses["C"] = 	0.012011;
        AtomicMasses["N"] = 	0.01400674;
        AtomicMasses["O"] = 	0.0159994;
        AtomicMasses["F"] = 	0.0189984032;
        AtomicMasses["Ne"] = 	0.0201797;
        AtomicMasses["Na"] = 	0.022989768;
        AtomicMasses["Mg"] = 	0.024305;
        AtomicMasses["Al"] = 	0.026981539;
        AtomicMasses["Si"] = 	0.0280855;
        AtomicMasses["P"] = 	0.030973762;
        AtomicMasses["S"] = 	0.032066;
        AtomicMasses["Cl"] = 	0.0354527;
        AtomicMasses["Ar"] = 	0.039948;
        AtomicMasses["K"] = 	0.0390983;
        AtomicMasses["Ca"] = 	0.040078;
        AtomicMasses["Sc"] = 	0.04495591;
        AtomicMasses["Ti"] = 	0.04788;
        AtomicMasses["V"] = 	0.0509415;
        AtomicMasses["Cr"] = 	0.0519961;
        AtomicMasses["Mn"] = 	0.05493805;
        AtomicMasses["Fe"] = 	0.055847;
        AtomicMasses["Co"] = 	0.0589332;
        AtomicMasses["Ni"] = 	0.0586934;
        AtomicMasses["Cu"] = 	0.063546;
        AtomicMasses["Zn"] = 	0.06539;
        AtomicMasses["Ga"] = 	0.069723;
        AtomicMasses["Ge"] = 	0.07261;
        AtomicMasses["As"] = 	0.07492159;
        AtomicMasses["Se"] = 	0.07896;
        AtomicMasses["Br"] = 	0.079904;
        AtomicMasses["Kr"] = 	0.0838;
        AtomicMasses["Rb"] = 	0.0854678;
        AtomicMasses["Sr"] = 	0.08762;
        AtomicMasses["Y"] = 	0.08890585;
        AtomicMasses["Zr"] = 	0.091224;
        AtomicMasses["Nb"] = 	0.09290638;
        AtomicMasses["Mo"] = 	0.09594;
        AtomicMasses["Tc"] = 	-0.0979072;
        AtomicMasses["Ru"] = 	0.10107;
        AtomicMasses["Rh"] = 	0.1029055;
        AtomicMasses["Pd"] = 	0.10642;
        AtomicMasses["Ag"] = 	0.1078682;
        AtomicMasses["Cd"] = 	0.112411;
        AtomicMasses["In"] = 	0.114818;
        AtomicMasses["Sn"] = 	0.11871;
        AtomicMasses["Sb"] = 	0.121757;
        AtomicMasses["Te"] = 	0.1276;
        AtomicMasses["I"] = 	0.12690447;
        AtomicMasses["Xe"] = 	0.13129;
        AtomicMasses["Cs"] = 	0.13290543;
        AtomicMasses["Ba"] = 	0.137327;
        AtomicMasses["La"] = 	0.1389055;
        AtomicMasses["Ce"] = 	0.140115;
        AtomicMasses["Pr"] = 	0.14090765;
        AtomicMasses["Nd"] = 	0.14424;
        AtomicMasses["Pm"] = 	-0.1449127;
        AtomicMasses["Sm"] = 	0.15036;
        AtomicMasses["Eu"] = 	0.151965;
        AtomicMasses["Gd"] = 	0.15725;
        AtomicMasses["Tb"] = 	0.15892534;
        AtomicMasses["Dy"] = 	0.1625;
        AtomicMasses["Ho"] = 	0.16493032;
        AtomicMasses["Er"] = 	0.16726;
        AtomicMasses["Tm"] = 	0.16893421;
        AtomicMasses["Yb"] = 	0.17304;
        AtomicMasses["Lu"] = 	0.174967;
        AtomicMasses["Hf"] = 	0.17849;
        AtomicMasses["Ta"] = 	0.1809479;
        AtomicMasses["W"] = 	0.18384;
        AtomicMasses["Re"] = 	0.186207;
        AtomicMasses["Os"] = 	0.19023;
        AtomicMasses["Ir"] = 	0.19222;
        AtomicMasses["Pt"] = 	0.19508;
        AtomicMasses["Au"] = 	0.19696654;
        AtomicMasses["Hg"] = 	0.20059;
        AtomicMasses["Tl"] = 	0.2043833;
        AtomicMasses["Pb"] = 	0.2072;
        AtomicMasses["Bi"] = 	0.20898037;
        AtomicMasses["Po"] = 	-0.2089824;
        AtomicMasses["At"] = 	-0.2099871;
        AtomicMasses["Rn"] = 	-0.2220176;
        AtomicMasses["Fr"] = 	-0.2230197;
        AtomicMasses["Ra"] = 	-0.2260254;
        AtomicMasses["Ac"] = 	-0.2270278;
        AtomicMasses["Th"] = 	0.2320381;
        AtomicMasses["Pa"] = 	0.23103588;
        AtomicMasses["U"] = 	0.2380289;
        AtomicMasses["Np"] = 	-0.2370482;
        AtomicMasses["Pu"] = 	-0.2440642;
        AtomicMasses["Am"] = 	-0.2430614;
        AtomicMasses["Cm"] = 	-0.2470703;
        AtomicMasses["Bk"] = 	-0.2470703;
        AtomicMasses["Cf"] = 	-0.2510796;
        AtomicMasses["Es"] = 	-0.252083;
        AtomicMasses["Fm"] = 	-0.2570951;
        AtomicMasses["Md"] = 	-0.2580984;
        AtomicMasses["No"] = 	-0.2591011;
        AtomicMasses["Lr"] = 	-0.2621098;
        AtomicMasses["Rf"] = 	-0.2611089;
        AtomicMasses["Ha"] = 	-0.2621144;
        AtomicMasses["Sg"] = 	-0.2631186;
        AtomicMasses["Ns"] = 	-0.2621231;
        AtomicMasses["Hs"] = 	-0.2651306;
        AtomicMasses["Mt"] = 	-0.2661378;
        AtomicMasses["Unn"] = 	-0.268;
        AtomicMasses["Unu"] = 	-0.269;


    }

};

#endif // DATA_H
