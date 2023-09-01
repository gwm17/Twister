#include "NuclearMap.h"

#include <fstream>

namespace Twister {

    NuclearMap::NuclearMap()
    {
        auto massFile = std::ifstream("etc/amdc_2020.txt");

        std::string junk;
        std::string elementSymbol;
        std::getline(massFile, junk);

        NucleusData data;
        uint32_t id;
        while (massFile >> data.Z)
        {
            massFile >> data.A;
            massFile >> elementSymbol;
            massFile >> data.atomicMassU;

            data.atomicMass = data.atomicMassU * s_u2MeV;
            data.isotopicMassU = data.atomicMassU - data.Z * s_eMass;
            data.isotopicMass = data.isotopicMassU * s_u2MeV;
            data.isotopicSymbol = std::to_string(data.A) + elementSymbol;

            id = GenerateID(data.Z, data.A);
            m_map[id] = data;
        }
    }

    NuclearMap::~NuclearMap()
    {
    }

    const NucleusData& NuclearMap::GetData(uint32_t Z, uint32_t A) const
    {
        uint32_t id = GenerateID(Z, A);
        auto iter = m_map.find(id);
        if (iter == m_map.end())
        {
            return s_invalidData;
        }
        return iter->second;
    }
}