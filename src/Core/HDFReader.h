#pragma once

#include "Cluster.h"
#include "highfive/H5File.hpp"

#include <filesystem>

namespace Twister {

    class HDFReader
    {
    public:
        HDFReader(const std::filesystem::path& path);
        ~HDFReader();

        const bool IsValid() const { return m_isValid; }

        std::pair<int, int> GetEventRange();
        int GetNumberOfClusters(int event);

        Cluster GetCluster(int event, int clusterIndex);

    private:
        HighFive::File m_hdfFile;
        HighFive::Group m_topGroup;
        bool m_isValid;
        
        static constexpr char s_clusterTopLevelGroup[] = "cluster";
        static constexpr char s_eventPrefix[] = "event";
        static constexpr char s_localClusterPrefix[] = "cluster";
    };
}