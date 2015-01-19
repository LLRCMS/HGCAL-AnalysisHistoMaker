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

#include "TVector2.h"

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


/*****************************************************************/
float SuperCluster::weightedEta() const
/*****************************************************************/
{
    //cout<<"Computing weighted Eta\n";
    double eta = 0.;
    for(const auto cluster : clusters())
    {
        //cout<<"  eta cluster = "<<cluster->eta()<<"\n";
        eta += (double)cluster->eta()*cluster->calibratedEt();
    }
    //cout<<" eta SC = "<<eta/m_et<<"\n";
    return eta/m_et;
}

/*****************************************************************/
float SuperCluster::weightedPhi() const
/*****************************************************************/
{
    //cout<<"Computing weighted Phi\n";
    double phi = 0.;
    double energySum = 0.;
    for(const auto cluster : clusters())
    {
        //cout<<"  phi cluster = "<<cluster->phi()<<"\n";
        double previousPhi = (energySum>0. ? phi/energySum : 0.);
        double diff = TVector2::Phi_mpi_pi(cluster->phi() - previousPhi);
        phi += (double)(previousPhi + diff)*cluster->calibratedEt();
        energySum += cluster->calibratedEt();
    }
    //cout<<" phi SC = "<<TVector2::Phi_mpi_pi(phi/m_et)<<"\n";
    return TVector2::Phi_mpi_pi(phi/m_et);
}
