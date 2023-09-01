#pragma once

#include <unordered_map>
#include <cstdint>
#include <string>

namespace Twister {

    struct NucleusData
    {
        uint32_t Z = 0;
        uint32_t A = 0;
        std::string isotopicSymbol = "Invalid";
        double isotopicMass = 0.0;
        double atomicMass = 0.0;
        double isotopicMassU = 0.0;
        double atomicMassU = 0.0;
    };

    class NuclearMap
    {
    public:
        NuclearMap();
        ~NuclearMap();

        const NucleusData& GetData(uint32_t Z, uint32_t A) const;

    private:
        //Using Szudzik pairing function to make unique id for two unsigned ints
        static constexpr uint32_t GenerateID(uint32_t Z, uint32_t A)
        {
    			return Z >= A ? (Z * Z + Z + A) : (A * A + Z);
        }
        std::unordered_map<uint32_t, NucleusData> m_map;

        //constants
        static constexpr double s_u2MeV = 931.4940954;
        static constexpr double s_eMass = 0.000548579909;
        NucleusData s_invalidData = NucleusData(); //Not literally const or static, but used functionally equivalently
    };
}