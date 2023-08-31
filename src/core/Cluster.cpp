#include "Cluster.h"

namespace Twister {

    Cluster::Cluster(const std::vector<std::vector<double>>& cloud, int label, int event, int index) :
        m_cloud(cloud), m_label(label), m_event(event), m_index(index)
    {
    }

    Cluster::~Cluster()
    {
    }

}