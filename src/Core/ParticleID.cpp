#include "ParticleID.h"

#include <iostream>
#include <fstream>

namespace Twister {

    ParticleID::ParticleID(const std::filesystem::path& path, const NuclearMap& nucMap) :
        m_isValid(false)
    {
        ReadFile(path, nucMap);
    }

    ParticleID::~ParticleID()
    {
    }

    void ParticleID::ReadFile(const std::filesystem::path& path, const NuclearMap& nucMap)
    {
        if (!std::filesystem::exists(path) || std::filesystem::is_directory(path))
        {
            std::cerr << "Filepath " << path << " given to ParticleID does not exist!" << std::endl;
            return;
        }

        auto jsonFile = std::ifstream(path);

        nlohmann::json jsonData = nlohmann::json::parse(jsonFile);

        uint32_t Z = jsonData["Z"];
        uint32_t A = jsonData["A"];
        m_particleData = nucMap.GetData(Z, A);
        if (m_particleData.isotopicSymbol == "Invalid")
        {
            std::cerr << "ParticleID found invalid nucleus with Z: " << Z << " A: " << A << std::endl;
            return;
        }

        std::vector<std::vector<double>> verticies = jsonData["verticies"];

        for (auto& vertex : verticies)
        {
            m_xPoints.push_back(vertex[0]);
            m_yPoints.push_back(vertex[1]);
        }

        m_name = jsonData["name"];
    }


    std::vector<bool> ParticleID::IsInside(const std::vector<Guess>& guesses) const
    {
        std::vector<bool> isInside;
        isInside.resize(guesses.size(), false);
        double slope;
        double x, y;
        for (std::size_t row=0; row < guesses.size(); row++)
        {
            slope = 0.0;
            x = guesses[row].dEdx;
            y = guesses[row].brho;
            for (std::size_t i=0; i<(m_xPoints.size() - 1); i++)
            {
                if (x == m_xPoints[i+1] && y == m_yPoints[i+1])
                {
                    isInside[row] = true;
                    break;
                }
                else if ((m_yPoints[i+1] > y) != (m_yPoints[i] > y))
                {
                    slope = (x - m_xPoints[i+1]) * (m_yPoints[i] - m_yPoints[i+1]) - (m_xPoints[i] - m_xPoints[i+1]) * (y - m_yPoints[i+1]);
                    if (slope == 0.0)
                    {
                        isInside[row] = true;
                        break;
                    }
                    else if ((slope < 0.0) != (m_yPoints[i] < m_yPoints[i+1]))
                    {
                        isInside[row] = !isInside[row];
                    }
                }
            }
        }

        return isInside;
    }
}