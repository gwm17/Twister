#include "Core/Solver.h"
#include "Core/ParticleID.h"
#include "Core/GuessReader.h"
#include "Core/NuclearMap.h"
#include "Core/Target.h"
#include "Utils/Timer.h"

#include <iostream>
#include <fstream>


int main(int argc, const char** argv)
{
    const std::filesystem::path guessFileName = "/media/gordon/ThesisData/NewData/Analyzed/a1975/estimates/run_0004.csv";
    const std::filesystem::path cloudFileName = "/media/gordon/ThesisData/NewData/Analyzed/a1975/clusters/run_0004.csv";
    const std::filesystem::path pidFileName = "/media/gordon/ThesisData/NewData/Analyzed/a1975/pid/ede_cut.json";

    Twister::NuclearMap nucMap;

    std::cout << "Loading particle id..." << std::endl;
    Twister::ParticleID pid(pidFileName, nucMap);

    std::cout << "Loading guesses..." << std::endl;
    auto validGuesses = Twister::ReadGuessesWithPID(guessFileName, pid);

    std::size_t row = 100;

    auto& chosenGuess = validGuesses[row];

    double magneticField = 2.9; //Tm
    double electricField = 60000.0; //Vm
    double density = 0.0000625; //Gas density

    std::cout << "Loading target..." << std::endl;
    Twister::Target target({1}, {2}, {2}, density, nucMap);

    std::cout << "Loading solver..." << std::endl;
    Twister::Solver solver(magneticField, electricField, target, pid.GetParticleData());

    std::cout << "Creating timesteps..." << std::endl;
    double t0 = 0.0;
    double tFinal= 1.0e-6;
    double tStep = 1.0e-9;
    int nSteps = (tFinal - t0) / tStep;
    std::vector<double> times;
    std::vector<std::vector<double>> results;
    results.resize(nSteps);
    for (int i=0; i<nSteps; i++)
    {
        times.push_back(t0 + (i+1) * tStep);
        results[i].resize(6, 0.);
    }

    std::cout << "Running solver..." << std::endl;
    Twister::Timer timer("Solver");
    solver.Run(chosenGuess, times, results);
    timer.Stop();

    std::cout << "Writing result..." << std::endl;
    std::ofstream output("test.csv");
    output << "x,y,z,vx,vy,vz" << std::endl;
    for(auto& point : results)
    {
        output << point[0] << "," << point[1] << "," << point[2] << "," << point[3] << "," << point[4] << "," << point[5] << std::endl;
    }

    std::cout << "Finished." << std::endl;
}