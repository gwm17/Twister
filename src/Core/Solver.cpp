#include "Solver.h"
#include "Constants.h"
#include "Precision.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

#include <iostream>

namespace Twister {

    int EquationsOfMotion(double t, const double y[], double dydt[], void* params)
    {
        (void)(t); //Ignore unused param warning
        EOMParams* paramHandle = static_cast<EOMParams*>(params);

        double speed = std::sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
        double unitVec[3] = {y[3]/speed, y[4]/speed, y[5]/speed};

        double kineticEnergy = paramHandle->ejectile.isotopicMass * (1.0/std::sqrt(1.0 - std::pow(speed / Constants::speedOfLight, 2.0)) - 1.0); //MeV
        //dE/dx -> MeV/g/cm^2 * J/MeV -> J/g/cm^2 * g/cm^3 -> J/cm * cm/m -> kg m/s^2 * 1/kg -> m/s^2
        double decel = paramHandle->target.GetAccelerationSI(paramHandle->ejectile, kineticEnergy);
        double qPerM = paramHandle->ejectile.Z * Constants::unitCharge / (paramHandle->ejectile.isotopicMass * Constants::MeV2kg); // C/kg

        //dr/dt = v
        dydt[0] = y[3]; //drx/dt
        dydt[1] = y[4]; //dry/dt
        dydt[2] = y[5]; //drz/dt
        //dv/dt = q/m * (E + vxB) - stopping
        dydt[3] = qPerM * paramHandle->Bfield * y[4] - decel * unitVec[0]; //dvx/dt
        dydt[4] = qPerM * -1.0 * paramHandle->Bfield * y[3] - decel * unitVec[1]; //dvy/dt
        dydt[5] = qPerM * paramHandle->Efield - decel * unitVec[2]; //dvz/dt

        return GSL_SUCCESS;
    }

    int Jacobian(double t, const double state[], double* dfdy, double dfdt[], void* params)
    {
        (void)(t); //Ignore unused param warning
        EOMParams* paramHandle = static_cast<EOMParams*>(params);
        gsl_matrix_view dfdyView = gsl_matrix_view_array(dfdy, 6, 6);
        gsl_matrix* dfdyMatrix = &dfdyView.matrix;

        double qPerM = paramHandle->ejectile.Z * Constants::unitCharge / (paramHandle->ejectile.isotopicMass * Constants::MeV2kg); // C/kg

        //df/dy

        //dfrx/dy
        gsl_matrix_set(dfdyMatrix, 0, 0, 0.0); //dfrx/drx
        gsl_matrix_set(dfdyMatrix, 0, 1, 0.0); //dfrx/dry
        gsl_matrix_set(dfdyMatrix, 0, 2, 0.0); //dfrx/drz
        gsl_matrix_set(dfdyMatrix, 0, 3, 1.0); //dfrx/dvx
        gsl_matrix_set(dfdyMatrix, 0, 4, 0.0); //dfrx/dvy
        gsl_matrix_set(dfdyMatrix, 0, 5, 0.0); //dfrx/dvz

        //dfry/dy
        gsl_matrix_set(dfdyMatrix, 1, 0, 0.0); //dfry/drx
        gsl_matrix_set(dfdyMatrix, 1, 1, 0.0); //dfry/dry
        gsl_matrix_set(dfdyMatrix, 1, 2, 0.0); //dfry/drz
        gsl_matrix_set(dfdyMatrix, 1, 3, 0.0); //dfry/dvx
        gsl_matrix_set(dfdyMatrix, 1, 4, 1.0); //dfry/dvy
        gsl_matrix_set(dfdyMatrix, 1, 5, 0.0); //dfry/dvz

        //dfrz/dy
        gsl_matrix_set(dfdyMatrix, 2, 0, 0.0); //dfrz/drx
        gsl_matrix_set(dfdyMatrix, 2, 1, 0.0); //dfrz/dry
        gsl_matrix_set(dfdyMatrix, 2, 2, 0.0); //dfrz/drz
        gsl_matrix_set(dfdyMatrix, 2, 3, 0.0); //dfrz/dvx
        gsl_matrix_set(dfdyMatrix, 2, 4, 0.0); //dfrz/dvy
        gsl_matrix_set(dfdyMatrix, 2, 5, 1.0); //dfrz/dvz

        //dfvx/dy
        gsl_matrix_set(dfdyMatrix, 3, 0, 0.0); //dfvx/drx
        gsl_matrix_set(dfdyMatrix, 3, 1, 0.0); //dfvx/dry
        gsl_matrix_set(dfdyMatrix, 3, 2, 0.0); //dfvx/drz
        gsl_matrix_set(dfdyMatrix, 3, 3, 0.0); //dfvx/dvx
        gsl_matrix_set(dfdyMatrix, 3, 4, qPerM * paramHandle->Bfield); //dfvx/dvy
        gsl_matrix_set(dfdyMatrix, 3, 5, 0.0); //dfvx/dvz

        //dfvy/dy
        gsl_matrix_set(dfdyMatrix, 4, 0, 0.0); //dfvy/drx
        gsl_matrix_set(dfdyMatrix, 4, 1, 0.0); //dfvy/dry
        gsl_matrix_set(dfdyMatrix, 4, 2, 0.0); //dfvy/drz
        gsl_matrix_set(dfdyMatrix, 4, 3, -1.0 * qPerM * paramHandle->Bfield); //dfvy/dvx
        gsl_matrix_set(dfdyMatrix, 4, 4, 0.0); //dfvy/dvy
        gsl_matrix_set(dfdyMatrix, 4, 5, 0.0); //dfvy/dvz

        //dfvz/dy
        gsl_matrix_set(dfdyMatrix, 5, 0, 0.0); //dfvz/drx
        gsl_matrix_set(dfdyMatrix, 5, 1, 0.0); //dfvz/dry
        gsl_matrix_set(dfdyMatrix, 5, 2, 0.0); //dfvz/drz
        gsl_matrix_set(dfdyMatrix, 5, 3, 0.0); //dfvz/dvx
        gsl_matrix_set(dfdyMatrix, 5, 4, 0.0); //dfvz/dvy
        gsl_matrix_set(dfdyMatrix, 5, 5, 0.0); //dfvz/dvz

        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 0.0;
        dfdt[4] = 0.0;
        dfdt[5] = 0.0;

        return GSL_SUCCESS;
    }

    Solver::Solver(double Bfield, double Efield, const Target& target, const NucleusData& ejectile) :
        m_driver(nullptr)
    {
        m_eomParams.Bfield = -1.0 * Bfield;
        m_eomParams.Efield = -1.0 * Efield;
        m_eomParams.target = target;
        m_eomParams.ejectile = ejectile;

        m_system = {EquationsOfMotion, Jacobian, 6, &m_eomParams};
        //Specify precision in terms of  x,y,z not vx,vy,vz
        m_driver = gsl_odeiv2_driver_alloc_y_new(&m_system, gsl_odeiv2_step_rkf45, 1.0e-9, 1.0e-6, 0.0);
    }

    Solver::~Solver()
    {
        if (m_driver != nullptr)
            gsl_odeiv2_driver_free(m_driver);
    }

    void Solver::Run(const Guess& initialGuess, const std::vector<double>& distanceSteps, std::vector<std::vector<double>>& results)
    {
        gsl_odeiv2_driver_reset(m_driver);

        if (distanceSteps.size() != results.size())
        {
            std::cerr << "Error at Solver::Run, the time vector and the result vector do not have the same length! Pre-allocate the result vector." << std::endl;
            return;
        }

        if (results[0].size() != s_sizeOfState)
        {
            std::cerr << "Error at Solver::Run, the result vector is not large enough to hold the state! Each result is of size " << s_sizeOfState << ". Pre-allocate the result vector." << std::endl;
            return;
        }

        //Create the initial state from the guessed values
        double tInitial = 0.0;
        double tFinal = 0.0;
        double y[6] = {0., 0., 0., 0., 0., 0.}; //x, y, z, vx, vy, vz
        double speed = (initialGuess.brho * 10.0 * 100.0 * Constants::QBrho2P * m_eomParams.ejectile.Z) / m_eomParams.ejectile.isotopicMass; //brho T*m -> kG*cm
	    speed *= Constants::speedOfLight; //Convert to m/s

        y[0] = initialGuess.vertexX * 0.001; //Convert to m
        y[1] = initialGuess.vertexY * 0.001;
        y[2] = initialGuess.vertexZ * 0.001;
        y[3] = speed * std::sin(initialGuess.polar) * std::cos(initialGuess.azimuthal);
        y[4] = speed * std::sin(initialGuess.polar) * std::sin(initialGuess.azimuthal);
        y[5] = speed * std::cos(initialGuess.polar);

        double kineticEnergy;

        for (std::size_t i=0; i<distanceSteps.size(); i++)
        {
            speed = std::sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
            tFinal = tInitial + distanceSteps[i] / speed;
            //Solve the ode
            gsl_odeiv2_driver_apply(m_driver, &tInitial, tFinal, y);
            //Save the result
            for (int j=0; j<s_sizeOfState; j++)
                results[i][j] = y[j];
            
            speed = std::sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
            kineticEnergy = m_eomParams.ejectile.isotopicMass * (1.0/std::sqrt(1.0 - std::pow(speed / Constants::speedOfLight, 2.0)) - 1.0); //MeV
            if (Precision::IsFloatAlmostEqual(kineticEnergy, 0.0, 2.0e-6))
            {
                break;
            }
        }

    }

}
