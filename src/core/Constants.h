#pragma once

namespace Twister {

    namespace Constants {
        static constexpr double speedOfLight = 299792458.0; // m/s
        static constexpr double MeV2Joule = 1.0 / 6.241509074e18 * 1.0e6; // J/MeV
        static constexpr double MeV2kg = 1.782661921e-36 * 1.0e6; //kg/MeV/c^2
        static constexpr double QBrho2P = 1.0e-9 * speedOfLight; //MeV/(Z*kG*cm)
    }
}