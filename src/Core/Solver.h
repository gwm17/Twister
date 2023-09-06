#pragma once

#include "Target.h"
#include "NuclearMap.h"
#include "GuessReader.h"
#include "Cluster.h"

namespace Twister {

    struct EOMParams
    {
        double Efield = 0.0; //Vm
        double Bfield = 0.0; //Tm
        Target target = Target();
        NucleusData ejectile = NucleusData();
    };

    struct ObjectiveParams
    {
        EOMParams eomParams;
        std::vector<double> steps;
        std::vector<std::vector<double>> trajectoryStorage;
        std::vector<std::vector<double>> data;
    };

    class Solver
    {
    public:
        Solver(double Bfield, double Efield, const Target& target, const NucleusData& ejectile);
        ~Solver();

        void SolveSystem(const Guess& initialGuess, const std::vector<double>& steps);
        Guess OptimizeSystem(const Guess& initialGuess, const std::vector<double>& steps, const Cluster& cluster);
        const std::vector<std::vector<double>>& GetTrajectory() const { return m_params.trajectoryStorage; }

    private:
        ObjectiveParams m_params;

        static constexpr std::size_t s_sizeOfState = 6; // x,y,z,vx,vy,vz
    };
}