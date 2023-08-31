#pragma once

#include <cmath>

namespace Twister {

    namespace Precision {

        static constexpr bool IsFloatAlmostEqual(double lhs, double rhs, double epsilon)
        {
            return std::fabs(lhs - rhs) < epsilon;
        }

        static constexpr bool IsFloatLessOrAlmostEqual(double lhs, double rhs, double epsilon)
        {
            return IsFloatAlmostEqual(lhs, rhs, epsilon) || lhs < rhs;
        }

        static constexpr bool IsFloatGreaterOrAlmostEqual(double lhs, double rhs, double epsilon)
        {
            return IsFloatAlmostEqual(lhs, rhs, epsilon) || lhs > rhs;
        }
    }
}