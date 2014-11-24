/**
 *  @file  SuperCluster.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    20/11/2014
 *
 *  @internal
 *     Created :  20/11/2014
 * Last update :  20/11/2014 17:37:33
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"

#include <iostream>


using namespace std;
using namespace AnHiMa;


/*****************************************************************/
bool AnHiMa::superClusterSort(const SuperCluster* sc1, const SuperCluster* sc2)
/*****************************************************************/
{
    return sc1->et() > sc2->et();
};

/*****************************************************************/
SuperCluster::SuperCluster():
    m_eta(0.),
    m_phi(0.),
    m_energy(0.),
    m_et(0.)
/*****************************************************************/
{
}

/*****************************************************************/
SuperCluster::SuperCluster(const SuperCluster& sc):
    m_eta(sc.eta()),
    m_phi(sc.phi()),
    m_energy(sc.energy()),
    m_et(sc.et())
/*****************************************************************/
{
    for(const auto cluster : sc.clusters())
    {
        m_clusters.push_back(cluster);
    }
}


/*****************************************************************/
SuperCluster::~SuperCluster()
/*****************************************************************/
{

}


/*****************************************************************/
void SuperCluster::addCluster(const Tower* cluster)
/*****************************************************************/
{
    //cerr<<"  Adding cluster ";
    //cerr<<"   eta="<<cluster->eta()<<",phi="<<cluster->phi()<<",Et="<<cluster->calibratedEt()<<"\n";
    if(m_clusters.size()==0) 
    {
        m_eta = cluster->eta();
        m_phi = cluster->phi();
    }
    m_clusters.push_back(cluster);
    m_energy += cluster->calibratedEnergy();
    m_et += cluster->calibratedEt();
}


/*****************************************************************/
const Tower* SuperCluster::cluster(unsigned i) const
/*****************************************************************/
{
    if(m_clusters.size()>i) return m_clusters[i];
    return 0;
}
