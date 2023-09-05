#include "Core/Solver.h"
#include "Core/ParticleID.h"
#include "Core/GuessReader.h"
#include "Core/NuclearMap.h"
#include "Core/Target.h"
#include "Core/Cluster.h"
#include "Core/HDFReader.h"
#include "Utils/Timer.h"

#include <iostream>
#include <fstream>


int main(int argc, const char** argv)
{
    const std::filesystem::path guessFileName = "/media/gordon/ThesisData/NewData/Analyzed/a1975/estimates/run_0004.csv";
    const std::filesystem::path cloudFileName = "/media/gordon/ThesisData/NewData/Analyzed/a1975/clusters/run_0004.h5";
    const std::filesystem::path pidFileName = "/media/gordon/ThesisData/NewData/Analyzed/a1975/pid/ede_cut.json";

    Twister::NuclearMap nucMap;

    std::cout << "Loading particle id..." << std::endl;
    Twister::ParticleID pid(pidFileName, nucMap);

    std::cout << "Loading guesses..." << std::endl;
    auto validGuesses = Twister::ReadGuessesWithPID(guessFileName, pid);

    std::size_t row = 666;

    auto& chosenGuess = validGuesses[row];

    std::cout << "Event: " << chosenGuess.event << " Cluster Index: " << chosenGuess.clusterIndex << std::endl;

    double magneticField = 2.9; //Tm
    double electricField = 6000.0; //Vm
    double density = 5.47e-6; //Gas density

    std::cout << "Loading target..." << std::endl;
    Twister::Target target({1}, {2}, {2}, density, nucMap);

    std::cout << "Loading solver..." << std::endl;
    Twister::Solver solver(magneticField, electricField, target, pid.GetParticleData());
    
    std::cout << "Loading cluster data..." << std::endl;

    Twister::HDFReader clusterReader(cloudFileName);
    Twister::Cluster cluster = clusterReader.GetCluster(chosenGuess.event, chosenGuess.clusterIndex);

    std::cout << "Calculating distance steps..." << std::endl;

    cluster.ConvertCloudToMeters();
    std::vector<double> steps = cluster.GetDistanceSteps(chosenGuess);
    std::vector<std::vector<double>> results;
    results.resize(steps.size());
    for (auto& r : results)
    {
        r.resize(6, 0.0);
    }

    std::cout << "Running solver..." << std::endl;
    Twister::Timer timer("Solver");
    solver.Run(chosenGuess, steps, results);
    timer.Stop();

    std::cout << "Writing result..." << std::endl;
    std::ofstream output("test.csv");
    output << "x,y,z,vx,vy,vz" << std::endl;
    for(auto& point : results)
    {
        if (point[0]  == 0.0 && point[1] == 0.0 && point[2] == 0.0)
            continue;
        output << point[0] << "," << point[1] << "," << point[2] << "," << point[3] << "," << point[4] << "," << point[5] << std::endl;
    }

    std::cout << "Finished." << std::endl;
}
