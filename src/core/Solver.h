#pragma once

#include "Target.h"
#include "NuclearMap.h"
#include "GuessReader.h"

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
        Solver(Guess& guess, double Bfield, double Efield, const Target& target, const NucleusData& ejectile);
        ~Solver();

        void Run(const std::vector<double>& times);

    private:
        std::vector<double> ConvertGuessToInitialValue();

        Guess m_initialGuess;
        EOMParams m_eomParams;
    };
}