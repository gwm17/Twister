#pragma once

#include <filesystem>
#include "Guess.h"
#include "ParticleID.h"

namespace Twister {

    std::vector<Guess> ReadGuesses(const std::filesystem::path& path);

    std::vector<Guess> ReadGuessesWithPID(const std::filesystem::path& path, const ParticleID& pid);
}