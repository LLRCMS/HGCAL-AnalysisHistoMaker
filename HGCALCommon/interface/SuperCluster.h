/**
 *  @file  SuperCluster.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    20/11/2014
 *
 *  @internal
 *     Created :  20/11/2014
 * Last update :  20/11/2014 17:23:21
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef SUPERCLUSTER_H
#define SUPERCLUSTER_H

#include <vector>

#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"

namespace AnHiMa
{
    class SuperCluster
    {
        public:
            SuperCluster();
            SuperCluster(const SuperCluster& sc);
            ~SuperCluster();

            void addCluster(const Tower* cluster);
            //
            int nClusters() const {return m_clusters.size();}
            const Tower* cluster(unsigned i) const;
            const std::vector<const Tower*>& clusters() const {return m_clusters;}
            //
            float eta()     const {return m_eta;}
            float phi()     const {return m_phi;}
            float weightedEta() const;
            float weightedPhi() const;
            float energy()  const {return m_energy;}
            float et()      const {return m_et;}
            //
            void setEnergy(double energy) {m_energy = energy;}
            void setEt(double et) {m_et = et;}

        private:
            double m_eta;
            double m_phi;
            double m_energy;
            double m_et;
            std::vector<const Tower*> m_clusters;
    };

    bool superClusterSort(const SuperCluster* sc1, const SuperCluster* sc2);

}


#endif
