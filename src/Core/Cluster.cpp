#include "Cluster.h"
#include <cmath>

namespace Twister {

    Cluster::Cluster(const std::vector<std::vector<double>>& cloud, int label, int event, int index) :
        m_cloud(cloud), m_label(label), m_event(event), m_index(index), m_isMeters(false)
    {
    }

    Cluster::~Cluster()
    {
    }

    void Cluster::ConvertCloudToMeters()
    {
        static constexpr double mm2Meters = 0.001;
        if (m_isMeters)
            return;

        for (std::vector<double>& point : m_cloud)
        {
            point[0] *= mm2Meters;
            point[1] *= mm2Meters;
            point[2] *= mm2Meters;
        }
        
        m_isMeters = true;
    }

    //Get only x,y,z data in meters
    std::vector<std::vector<double>> Cluster::GetCloudPoints() const 
    {
        std::vector<std::vector<double>> points;
        points.resize(m_cloud.size());

        double scale = 1.0;
        if (!m_isMeters)
            scale = 0.001;

        for (std::size_t i=0; i<m_cloud.size(); i++)
        {
            points[i].resize(3);
            points[i][0] = m_cloud[i][0] * scale;
            points[i][1] = m_cloud[i][1] * scale;
            points[i][2] = m_cloud[i][2] * scale;
        }

        return points;
    }

    std::vector<double> Cluster::GetDistanceSteps(const Guess& guess) const
    {
        std::vector<double> steps;
        steps.resize(m_cloud.size(), 0.0);
        double dist=0.0;

        double scaleFactor = 1.0;
        if (m_isMeters)
        {
            scaleFactor = 0.001;
        }
        for (std::size_t i=0; i<m_cloud.size(); i++)
        {
            auto& nextPoint = m_cloud[i];
            if (i == 0)
            {
                dist = std::sqrt(std::pow((nextPoint[0] - guess.vertexX*scaleFactor), 2.0) + std::pow((nextPoint[1] - guess.vertexY*scaleFactor), 2.0) + std::pow((nextPoint[2] - guess.vertexZ*scaleFactor), 2.0));
            }
            else
            {
                auto& prevPoint = m_cloud[i-1];
                dist = std::sqrt(std::pow((nextPoint[0] - prevPoint[0]), 2.0) + std::pow((nextPoint[1] - prevPoint[1]), 2.0) + std::pow((nextPoint[2] - prevPoint[2]), 2.0));
            }
            steps[i] = dist;
        }
        return steps;
    }

}