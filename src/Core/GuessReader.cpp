#include "GuessReader.h"

#include <iostream>
#include <fstream>

#include "csv.hpp"

#define NUMBER_OF_COLUMNS = 16

namespace Twister {

    std::vector<Guess> ReadGuesses(const std::filesystem::path& path)
    {
        std::vector<Guess> result;
        if (!std::filesystem::exists(path) || std::filesystem::is_directory(path))
        {
            std::cerr << "File path " << path << " given to ReadGuesses does not exist!" << std::endl;
            return result;
        }

        Guess guess;
        csv::CSVReader reader(path.string());
        for (csv::CSVRow& row: reader)
        {
            guess.event = row["event"].get<int>();
            guess.clusterIndex = row["cluster_index"].get<int>();
            guess.clusterLabel = row["cluster_label"].get<int>();
            guess.vertexX = row["vertex_x"].get<double>();
            guess.vertexY = row["vertex_y"].get<double>();
            guess.vertexZ = row["vertex_z"].get<double>();
            guess.centerX = row["center_x"].get<double>();
            guess.centerY = row["center_y"].get<double>();
            guess.centerZ = row["center_z"].get<double>();
            guess.polar = row["polar"].get<double>();
            guess.azimuthal = row["azimuthal"].get<double>();
            guess.brho = row["brho"].get<double>();
            guess.dEdx = row["dEdx"].get<double>();
            guess.dE = row["dE"].get<double>();
            guess.arcLength = row["arclength"].get<double>();
            guess.direction = row["direction"].get<int>();

            result.push_back(guess);
        }
        return result;
    }

    std::vector<Guess> ReadGuessesWithPID(const std::filesystem::path& path, const ParticleID& pid)
    {
        std::vector<Guess> guesses = ReadGuesses(path);

        std::vector<Guess> validGuesses;
        std::vector<bool> isInside = pid.IsInside(guesses);
        for (std::size_t i=0; i<guesses.size(); i++)
        {
            if (isInside[i])
            {
                validGuesses.push_back(guesses[i]);
            }
        }

        return validGuesses;
    }
}