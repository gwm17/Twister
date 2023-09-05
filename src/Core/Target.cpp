#include "Target.h"
#include "Precision.h"
#include "Constants.h"

#include <iostream>

namespace Twister {

    //z,a: istope list of material compound, stoich: compound stoichometry, density: material density in g/cm^3
	Target::Target(const std::vector<uint32_t>& z, const std::vector<uint32_t>& a, const std::vector<int>& stoich, double density, NuclearMap& nucMap)
	{
		for(size_t i=0; i<z.size(); i++)
		{
			m_material.add_element(nucMap.GetData(z[i], a[i]).atomicMassU, z[i], stoich[i]);
		}

		m_material.density(density); //g/cm^3
	}

    Target::~Target()
    {
    }

    //Energy in MeV, returns MeV/g/cm^2
    double Target::dEdx(const NucleusData& projectile,  double energy)
    {
        if (Precision::IsFloatAlmostEqual(energy, 0.0, s_epsilon))
            return 0.0;

        catima::Projectile particle(projectile.isotopicMassU, projectile.Z);
        particle.T = energy / projectile.isotopicMassU;
        return catima::dedx(particle, m_material);
    }

    //Energy in MeV, returns accel in m/s^2
    double Target::GetAccelerationSI(const NucleusData& projectile, double energy)
    {
        double dedx = dEdx(projectile, energy);
        double massKg = projectile.isotopicMass * Constants::MeV2kg;
        //dE/dx -> MeV/g/cm^2 * J/MeV -> J/g/cm^2 * g/cm^3 -> J/cm * cm/m -> kg m/s^2 * 1/kg -> m/s^2
        return (dedx * Constants::MeV2Joule * GetDensity() * 100.0) / massKg;
    }

    double Target::GetPathLength(const NucleusData& projectile, double startEnergy, double finalEnergy)
    {
        double densityInv = 1.0/m_material.density();
		catima::Projectile proj(projectile.isotopicMassU, projectile.Z);
		proj.T = startEnergy / proj.A;
		double stopRange = catima::range(proj, m_material); //get the total range for startEnergy -> 0, returns g/cm^2!
		proj.T = finalEnergy/proj.A;
		double finalRange = catima::range(proj, m_material); //get the range bound from the final energy;
		return (stopRange - finalRange) * densityInv * 0.01; //convert to m
    }

    double Target::GetAngularStraggling(const NucleusData& projectile, double energy)
    {
        catima::Projectile proj(projectile.isotopicMassU, projectile.Z);
		proj.T = energy/proj.A;
		return catima::angular_straggling(proj, m_material);
    }

}
