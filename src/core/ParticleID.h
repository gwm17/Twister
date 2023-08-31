#pragma once

#include <filesystem>
#include <vector>
#include "GuessReader.h"

namespace Twister {

    class ParticleID
    {
    public:
        ParticleID(const std::filesystem::path& path);
        ~ParticleID();

        std::vector<bool> IsInside(const std::vector<Guess>& guesses) const;

    private:
        std::vector<double> m_xPoints;
        std::vector<double> m_yPoints;

    };
}