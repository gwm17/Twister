#include "Solver.h"
#include "Constants.h"
#include "Precision.h"

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>

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

    void Solve(const gsl_vector* initialGuess, const std::vector<double>& distanceSteps, std::vector<std::vector<double>>& trajectory, EOMParams& params)
    {
        static constexpr int sizeOfState = 6;

        gsl_odeiv2_system system  = {EquationsOfMotion, Jacobian, 6, &params};
        gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1.0e-9, 1.0e-6, 0.0);

        //Create the initial state from the guessed values
        double tInitial = 0.0;
        double tFinal = 0.0;
        double y[6] = {0., 0., 0., 0., 0., 0.}; //x, y, z, vx, vy, vz
        double polar = gsl_vector_get(initialGuess, 3);
        double azimuthal = gsl_vector_get(initialGuess, 4);
        double brho = gsl_vector_get(initialGuess, 5);
        double speed = (brho * 10.0 * 100.0 * Constants::QBrho2P * params.ejectile.Z) / params.ejectile.isotopicMass; //brho T*m -> kG*cm
	    speed *= Constants::speedOfLight; //Convert to m/s

        y[0] = gsl_vector_get(initialGuess, 0) * 0.001; //Convert to m
        y[1] = gsl_vector_get(initialGuess, 1) * 0.001;
        y[2] = gsl_vector_get(initialGuess, 2) * 0.001;
        y[3] = speed * std::sin(polar) * std::cos(azimuthal);
        y[4] = speed * std::sin(polar) * std::sin(azimuthal);
        y[5] = speed * std::cos(polar);

        double kineticEnergy;

        for (std::size_t i=0; i<distanceSteps.size(); i++)
        {
            speed = std::sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
            tFinal = tInitial + distanceSteps[i] / speed;
            //Solve the ode
            gsl_odeiv2_driver_apply(driver, &tInitial, tFinal, y);
            //Save the result
            for (int j=0; j<sizeOfState; j++)
                trajectory[i][j] = y[j];
            
            speed = std::sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]);
            kineticEnergy = params.ejectile.isotopicMass * (1.0/std::sqrt(1.0 - std::pow(speed / Constants::speedOfLight, 2.0)) - 1.0); //MeV
            if (Precision::IsFloatAlmostEqual(kineticEnergy, 0.0, 2.0e-6))
            {
                break;
            }
        }

        gsl_odeiv2_driver_free(driver);
    }

    double Objective(const gsl_vector* guess, void* params)
    {
        ObjectiveParams* oParams = static_cast<ObjectiveParams*>(params);
        Solve(guess, oParams->steps, oParams->trajectoryStorage, oParams->eomParams);
        double sumResiduals = 0.0;
        double n = oParams->trajectoryStorage.size();
        for (std::size_t i=0; i<oParams->trajectoryStorage.size(); i++)
        {
            auto& trajPoint = oParams->trajectoryStorage[i];
            auto& dataPoint = oParams->data[i];
            sumResiduals += std::sqrt(std::pow(trajPoint[0] - dataPoint[0], 2.0) + std::pow(trajPoint[1] - dataPoint[1], 2.0) + std::pow(trajPoint[2] - dataPoint[2], 2.0));
        }
        return sumResiduals / n;
    }

    std::size_t Minimize(gsl_vector* guess, const ObjectiveParams& params)
    {
        static constexpr int sizeOfParameters = 6;
        static constexpr std::size_t maxIters = 1000;
        const gsl_multimin_fminimizer_type* type = gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_function function;

        function.n = sizeOfParameters;
        function.f = Objective;
        function.params = (void*) &params; // gross

        gsl_vector* parameterStepSizes = gsl_vector_alloc(sizeOfParameters);
        gsl_vector_set(parameterStepSizes, 0, 0.001); //1mm
        gsl_vector_set(parameterStepSizes, 1, 0.001); //1mm
        gsl_vector_set(parameterStepSizes, 2, 0.001); //1mm
        gsl_vector_set(parameterStepSizes, 3, 0.09); //~5degrees
        gsl_vector_set(parameterStepSizes, 4, 0.09); //~5degrees
        gsl_vector_set(parameterStepSizes, 5, 0.2); //Tm

        gsl_multimin_fminimizer* minimizer = gsl_multimin_fminimizer_alloc(type, sizeOfParameters);
        gsl_multimin_fminimizer_set(minimizer, &function, guess, parameterStepSizes);

        double size;
        int status = GSL_CONTINUE;
        std::size_t nIter = 0;
        while(status == GSL_CONTINUE && nIter < maxIters)
        {
            nIter++;
            //Iterate
            status = gsl_multimin_fminimizer_iterate(minimizer);

            if (status)
                break;

            size = gsl_multimin_fminimizer_size(minimizer);
            status = gsl_multimin_test_size(size, 1.0e-2);
        }

        gsl_vector_free(parameterStepSizes);
        gsl_multimin_fminimizer_free(minimizer);

        if (status == GSL_SUCCESS)
        {
            std::cout << "Minimizer successfully converged!" << std::endl;
        }
        else
        {
            std::cerr << "Minimizer did not converge! GSL code: " << status << std::endl;
        }

        return nIter;
    }

    Solver::Solver(double Bfield, double Efield, const Target& target, const NucleusData& ejectile)
    {
        m_params.eomParams.Bfield = -1.0 * Bfield;
        m_params.eomParams.Efield = -1.0 * Efield;
        m_params.eomParams.target = target;
        m_params.eomParams.ejectile = ejectile;
    }

    Solver::~Solver()
    {
    }

    void Solver::SolveSystem(const Guess& initialGuess, const std::vector<double>& steps)
    {
        gsl_vector* guessVector = gsl_vector_alloc(s_sizeOfState);
        gsl_vector_set(guessVector, 0, initialGuess.vertexX * 0.001);
        gsl_vector_set(guessVector, 1, initialGuess.vertexY * 0.001);
        gsl_vector_set(guessVector, 2, initialGuess.vertexZ * 0.001);
        gsl_vector_set(guessVector, 3, initialGuess.polar);
        gsl_vector_set(guessVector, 4, initialGuess.azimuthal);
        gsl_vector_set(guessVector, 5, initialGuess.brho);

        m_params.steps = steps;
        m_params.trajectoryStorage.resize(steps.size());
        for (auto& point : m_params.trajectoryStorage)
            point.resize(s_sizeOfState, 0.0);

        Solve(guessVector, m_params.steps, m_params.trajectoryStorage, m_params.eomParams);
    }

    Guess Solver::OptimizeSystem(const Guess& initialGuess, const std::vector<double>& steps, const Cluster& data)
    {
        gsl_vector* guessVector = gsl_vector_alloc(s_sizeOfState);
        gsl_vector_set(guessVector, 0, initialGuess.vertexX * 0.001);
        gsl_vector_set(guessVector, 1, initialGuess.vertexY * 0.001);
        gsl_vector_set(guessVector, 2, initialGuess.vertexZ * 0.001);
        gsl_vector_set(guessVector, 3, initialGuess.polar);
        gsl_vector_set(guessVector, 4, initialGuess.azimuthal);
        gsl_vector_set(guessVector, 5, initialGuess.brho);
        m_params.steps = steps;
        m_params.data = data.GetCloudPoints();
        m_params.trajectoryStorage.resize(steps.size());
        for (auto& point : m_params.trajectoryStorage)
            point.resize(s_sizeOfState, 0.0);
        
        std::size_t nIter = Minimize(guessVector, m_params);

        std::cout << "Minimizer took " << nIter << " iterations." << std::endl;

        Guess bestGuess;
        bestGuess.vertexX = gsl_vector_get(guessVector, 0) * 1000.0;
        bestGuess.vertexY = gsl_vector_get(guessVector, 1) * 1000.0;
        bestGuess.vertexZ = gsl_vector_get(guessVector, 2) * 1000.0;
        bestGuess.polar = gsl_vector_get(guessVector, 3);
        bestGuess.azimuthal = gsl_vector_get(guessVector, 4);
        bestGuess.azimuthal = gsl_vector_get(guessVector, 5);

        gsl_vector_free(guessVector);
        return bestGuess;
    }

}
