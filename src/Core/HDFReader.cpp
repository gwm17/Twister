#include "HDFReader.h"

#include <iostream>
#include <string>

namespace Twister {

    HDFReader::HDFReader(const std::filesystem::path& path) :
        m_hdfFile(path.string(), HighFive::File::ReadOnly)
    {
        if (m_hdfFile.exist(s_clusterTopLevelGroup))
        {
            m_isValid = true;
            m_topGroup = m_hdfFile.getGroup(s_clusterTopLevelGroup);
        }
        else
            m_isValid = false;
    }

    HDFReader::~HDFReader()
    {
    }

    std::pair<int, int> HDFReader::GetEventRange()
    {
        int minEvent = m_topGroup.getAttribute("min_event").read<int>();
        int maxEvent = m_topGroup.getAttribute("max_event").read<int>();
        return std::make_pair(minEvent, maxEvent);
    }

    int HDFReader::GetNumberOfClusters(int event)
    {
        return m_topGroup.getGroup(s_eventPrefix + std::to_string(event)).getAttribute("nclusters").read<int>();
    }

    Cluster HDFReader::GetCluster(int event, int clusterIndex)
    {
        HighFive::Group localCluster = m_topGroup.getGroup(s_eventPrefix + std::string("_") + std::to_string(event)).getGroup(s_localClusterPrefix + std::string("_") + std::to_string(clusterIndex));
        int clusterLabel = localCluster.getAttribute("label").read<int>();
        std::vector<std::vector<double>> cloud = localCluster.getDataSet("cloud").read<std::vector<std::vector<double>>>();
        return Cluster(cloud, clusterLabel, event, clusterIndex);
    }
}