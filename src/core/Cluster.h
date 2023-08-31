#pragma once

#include <vector>
#include <array>

namespace Twister {

    class Cluster
    {
    public:
        Cluster(const std::vector<std::array<double, 6>>& cloud, int label, int event, int index);
        ~Cluster();

        const std::vector<std::array<double, 6>>& GetCloud() const { return m_cloud; }
        const int GetLabel() const { return m_label; }
        const int GetEvent() const { return m_event; }
        const int GetIndex() const { return m_index; }

    private:
        std::vector<std::array<double, 6>> m_cloud;
        int m_label;
        int m_event;
        int m_index;
    };
}