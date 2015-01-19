/**
 *  @file  TriggerJetAlgorithm.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    14/01/2015
 *
 *  @internal
 *     Created :  14/01/2015
 * Last update :  14/01/2015 13:59:18
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#include <fstream>
#include <list>

#include "TEnv.h"
#include "TFile.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerJetAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"

using namespace std;
using namespace AnHiMa;


/*****************************************************************/
TriggerJetAlgorithm::TriggerJetAlgorithm()
/*****************************************************************/
{
}



/*****************************************************************/
TriggerJetAlgorithm::~TriggerJetAlgorithm()
/*****************************************************************/
{
}



/*****************************************************************/
void TriggerJetAlgorithm::initialize(const EventHGCAL& event, TEnv& params)
/*****************************************************************/
{
}



/*****************************************************************/
void TriggerJetAlgorithm::hcalSeeding(const EventHGCAL& event, const vector<SimHit>& hits, vector<Tower>& seeds, double clusterSize)
/*****************************************************************/
{
    TowerCalibrator towerCalibrator;

    // find hits above 5 MIPs in FH
    list<const SimHit*> selectedHits;
    for(const auto& hit : hits)
    {
        // should be changed with PU dependent threshold
        if(hit.energy()>10.*mip_HF) selectedHits.push_back(&hit);
    }
    selectedHits.sort(AnHiMa::simHitSort);

    // Build clusters around selected hits, using only selected hits
    // If a hit is clustered, remove it from the list for the next clusters
    for(auto it1=selectedHits.begin(); it1!=std::prev(selectedHits.end()) && it1!=selectedHits.end(); it1++)
    {
        Tower cluster;
        const SimHit* hit1 = *it1;
        //cout<<"Main hit E="<<hit1->energy()<<", eta="<<hit1->eta()<<",phi="<<hit1->phi()<<"\n";
        cluster.addHit(*hit1);
        //
        GlobalPoint point1(hit1->x(), hit1->y(), hit1->z());
        GlobalPoint proj1( GlobalPoint::Polar(point1.theta(), point1.phi(), point1.mag()*354.77/point1.z()) );
        //cout<<"   Projected on eta="<<proj1.eta()<<",phi="<<proj1.phi()<<", z="<<proj1.z()<<"\n";
        for(auto it2=std::next(it1); it2!=selectedHits.end();)
        {
            const SimHit* hit2 = *it2;
            GlobalPoint point2(hit2->x(), hit2->y(), hit2->z());
            GlobalPoint proj2( GlobalPoint::Polar(point2.theta(), point2.phi(), point2.mag()*354.77/point2.z()) );
            double dx = proj1.x()-proj2.x();
            double dy = proj1.y()-proj2.y();
            double dr = sqrt(dx*dx+dy*dy);
            if(dr<clusterSize) // cluster hit
            {
                //cout<<"     Adding hit E="<<hit2->energy()<<", eta="<<hit2->eta()<<",phi="<<hit2->phi()<<"\n";
                cluster.addHit(*hit2);
                it2 = selectedHits.erase(it2);
            }
            else // move forward
            {
                it2++;
            }
        }
        // calibrate energy
        towerCalibrator.calibrate(cluster);
        //cout<<"Cluster energy="<<cluster.energy()<<", eta="<<cluster.eta()<<",phi="<<cluster.phi()<<"\n";
        // filter cluster candidates
        //if(cluster.calibratedEt()<1.) continue;
        // fill cluster
        seeds.push_back(cluster);
    }

}


/*****************************************************************/
void TriggerJetAlgorithm::superClustering(const EventHGCAL& event, const vector<Tower>& seeds, vector<SuperCluster>& superClusters, double coneSize)
/*****************************************************************/
{
    vector<SuperCluster*> sortedSuperClusters;
    // Find seed clusters
    // Seed clusters are clusters with max Et in a <coneSize> eta x phi cone
    for(const auto& cluster : seeds)
    {
        double eta0 = cluster.eta();
        double phi0 = cluster.phi();
        double et0 = cluster.et();
        bool isSeedCluster = true;
        // find neighbour clusters
        for(const auto& subcluster : seeds)
        {
            double eta = subcluster.eta();
            double phi = subcluster.phi();
            double et = subcluster.et();
            double deta = (eta-eta0);
            double dphi = TVector2::Phi_mpi_pi(phi-phi0);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(fabs(deta)<1.e-5 && fabs(dphi)<1.e-5) continue;
            if(dr>coneSize) continue;
            if(et>et0) 
            {
                isSeedCluster = false;
                break;
            }
        }
        if(isSeedCluster) 
        {
            superClusters.push_back( SuperCluster() );
            superClusters.back().addCluster(&cluster);
        }
    }
    for(auto& sc : superClusters) sortedSuperClusters.push_back(&sc);
    sort(sortedSuperClusters.begin(), sortedSuperClusters.end(), AnHiMa::superClusterSort);

    std::set<const Tower*> usedClusters;
    // Gather subclusters
    for(auto sc : sortedSuperClusters)
    {
        double eta0 = sc->eta();
        double phi0 = sc->phi();
        for(const auto& subcluster : seeds)
        {
            if(usedClusters.find(&subcluster)!=usedClusters.end()) continue;
            double eta = subcluster.eta();
            double phi = subcluster.phi();
            double deta = (eta-eta0);
            double dphi = TVector2::Phi_mpi_pi(phi-phi0);
            double dr = sqrt(deta*deta + dphi*dphi);
            // discard seed cluster
            if(fabs(deta)<1.e-5 && fabs(dphi)<1.e-5) continue;
            // outside window
            if(dr>coneSize) continue;
            //
            sc->addCluster(&subcluster);
            usedClusters.insert(&subcluster);
        }
    }
    //
}

/*****************************************************************/
void TriggerJetAlgorithm::coneClustering(const EventHGCAL& event, const vector<SimHit>& hits, const vector<SuperCluster>& superClusters, vector<Tower>& cones, double coneSize)
/*****************************************************************/
{
    TowerCalibrator towerCalibrator;

    for(const auto& sc : superClusters)
    {
        double eta0 = sc.eta();
        double phi0 = sc.phi();
        // first find barycenter of the ROI
        Tower cluster;
        for(const auto& hit: hits)
        {
            float eta = hit.eta();
            float phi = hit.phi();
            double deta = (eta-eta0);
            double dphi =  TVector2::Phi_mpi_pi(phi-phi0);
            float dr = sqrt(deta*deta+dphi*dphi);
            if(dr>0.4) continue;
            cluster.addHit(hit);
        }
        eta0 = cluster.eta();
        phi0 = cluster.phi();
        // then build cone around this barycenter 
        //Tower cone;
        //for(const auto& hit: hits)
        //{
        //    float eta = hit.eta();
        //    float phi = hit.phi();
        //    double deta = (eta-eta0);
        //    double dphi =  TVector2::Phi_mpi_pi(phi-phi0);
        //    float dr = sqrt(deta*deta+dphi*dphi);
        //    if(dr>coneSize) continue;
        //    cone.addHit(hit);
        //}
        towerCalibrator.calibrate(cluster);
        cones.push_back(cluster);
    }

}



