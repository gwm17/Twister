#pragma once

#include <filesystem>

namespace Twister {

    enum Direction
    {
        NONE=-1,
        FORWARD=0,
        BACKWARD=1
    };

    struct Guess
    {
        int event;
        int clusterIndex;
        int clusterLabel;
        int direction;
        double polar;
        double azimuthal;
        double brho;
        double dEdx;
        double dE;
        double arcLength;
        double centerX;
        double centerY;
        double centerZ;
        double vertexX;
        double vertexY;
        double vertexZ;
    };

    std::vector<Guess> ReadGuesses(const std::filesystem::path& path);
}