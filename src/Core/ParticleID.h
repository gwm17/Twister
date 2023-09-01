#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include "Guess.h"
#include "NuclearMap.h"
#include "nlohmann/json.hpp"

namespace Twister {
    
    class ParticleID
    {
    public:
        ParticleID(const std::filesystem::path& path, const NuclearMap& nucMap);
        ~ParticleID();

        std::vector<bool> IsInside(const std::vector<Guess>& guesses) const;
        const bool IsValid() const { return m_isValid; }
        const std::string& GetName() const { return m_name; }
        const NucleusData& GetParticleData() const { return m_particleData; }

    private:
        void ReadFile(const std::filesystem::path& path, const NuclearMap& nucMap);
        std::vector<double> m_xPoints;
        std::vector<double> m_yPoints;
        bool m_isValid;
        std::string m_name;

        NucleusData m_particleData;
    };
}