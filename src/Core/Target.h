#pragma once

#include <string>
#include <vector>
#include <cmath>
#include "catima/catima.h"
#include "NuclearMap.h"

namespace Twister {

	class Target {
	
	public:
        Target() = default;
	 	Target(const std::vector<uint32_t>& z, const std::vector<uint32_t>& a, const std::vector<int>& stoich, double density, NuclearMap& nucMap);
	 	~Target();

	 	double dEdx(const NucleusData& data, double energy);
		//Takes in energy in MeV -> returns m/s^2
		double GetAccelerationSI(const NucleusData& data, double energy);
		double GetPathLength(const NucleusData& data, double startEnergy, double finalEnergy); //Returns pathlength for a particle w/ startE to reach finalE (cm)
		double GetAngularStraggling(const NucleusData& data, double energy); //Returns planar angular straggling in radians for a particle with energy
	 	inline double GetDensity() { return m_material.density(); } //g/cm^3
	
	private:
		catima::Material m_material;

		static constexpr double s_epsilon = 1.0e-6;
	};

}