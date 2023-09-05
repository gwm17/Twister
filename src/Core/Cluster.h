#pragma once

#include <vector>
#include <array>
#include "Guess.h"

namespace Twister {

    class Cluster
    {
    public:
        Cluster(const std::vector<std::vector<double>>& cloud, int label, int event, int index);
        ~Cluster();

        const std::vector<std::vector<double>>& GetCloud() const { return m_cloud; }
        const int GetLabel() const { return m_label; }
        const int GetEvent() const { return m_event; }
        const int GetIndex() const { return m_index; }
        std::vector<double> GetDistanceSteps(const Guess& guess) const;
        void ConvertCloudToMeters();


    private:
        std::vector<std::vector<double>> m_cloud;
        int m_label;
        int m_event;
        int m_index;

        bool m_isMeters;

        static constexpr std::size_t s_pointSize = 6; // Length of the second axis (length of a single point)
    };
}