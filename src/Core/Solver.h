#pragma once

#include "Target.h"
#include "NuclearMap.h"
#include "GuessReader.h"

#include <gsl/gsl_odeiv2.h>

namespace Twister {

    struct EOMParams
    {
        double Efield = 0.0; //Vm
        double Bfield = 0.0; //Tm
        Target target = Target();
        NucleusData ejectile = NucleusData();
    };

    class Solver
    {
    public:
        Solver(double Bfield, double Efield, const Target& target, const NucleusData& ejectile);
        ~Solver();

        void Run(const Guess& initialGuess, const std::vector<double>& times, std::vector<std::vector<double>>& results);

    private:
        std::vector<double> ConvertGuessToInitialValue();

        EOMParams m_eomParams;

        gsl_odeiv2_system m_system;
        gsl_odeiv2_driver* m_driver;

        static constexpr std::size_t s_sizeOfState = 6; // x,y,z,vx,vy,vz
    };
}