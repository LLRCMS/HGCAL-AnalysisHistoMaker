/**
 *  @file  AnalysisClustering.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 13:44:32
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/StudyClustering/interface/AnalysisClustering.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisClustering::AnalysisClustering():IAnalysis(),
    m_genParticle(0),
    m_matchedSuperCluster(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisClustering::~AnalysisClustering()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisClustering::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // Read parameters
    //string pileupParamsFile = m_reader.params().GetValue("PileupParams", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/ElectronClusterThreshold/data/thresholdParameters.txt");

    m_egammaAlgo.initialize(event(), m_reader.params());
    return true;
}


/*****************************************************************/
void AnalysisClustering::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_clusters.clear();
    m_clustersHard.clear();
    m_clustersHardNoPUThreshold.clear();
    m_superClusters.clear();
    m_superClustersHard.clear();
    m_superClustersHardNoPUThreshold.clear();
    //cerr<<"seeding\n";
    m_egammaAlgo.seeding(event(), m_seeds);
    //cerr<<"clustering\n";
    m_egammaAlgo.clustering(event(), m_seeds, m_clusters);
    m_egammaAlgo.clustering(event(), m_seeds, m_clustersHard, true, true);
    m_egammaAlgo.clustering(event(), m_seeds, m_clustersHardNoPUThreshold, true, false);
    // calibration
    m_egammaAlgo.clusterCorrection(m_clusters);
    //cerr<<"superClustering\n";
    m_egammaAlgo.superClustering(m_clusters, m_superClusters);
    m_egammaAlgo.superClustering(m_clustersHard, m_superClustersHard);
    m_egammaAlgo.superClustering(m_clustersHardNoPUThreshold, m_superClustersHardNoPUThreshold);


    m_genParticle = &event().genparticles()[0];
    if(fabs(m_genParticle->Eta())<2.9 && fabs(m_genParticle->Eta())>1.6) 
    {
        fillHistos();
    }
}

/*****************************************************************/
void AnalysisClustering::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;



    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //
    m_histos.FillHisto(10+hoffset, m_seeds.size(), weight, sysNum);

    for(const auto& tower : m_seeds)
    {
        m_histos.FillHisto(11+hoffset, fabs(tower.eta()), weight, sysNum);
        m_histos.FillHisto(12+hoffset, tower.phi(), weight, sysNum);
        m_histos.FillHisto(13+hoffset, fabs(tower.eta()), tower.phi(), weight, sysNum);
        m_histos.FillHisto(14+hoffset, tower.calibratedEt(), weight, sysNum);

        //m_histos.Fill1BinHisto(100+hoffset, event().event()-1, fabs(tower.eta()), tower.phi(), weight, sysNum);
    }

    for(const auto& cluster : m_clusters)
    {
        m_histos.FillHisto(21+hoffset, fabs(cluster.eta()), weight, sysNum);
        m_histos.FillHisto(22+hoffset, cluster.phi(), weight, sysNum);
        m_histos.FillHisto(24+hoffset, cluster.calibratedEt(), weight, sysNum);
    }

    for(int i=0;i<2;i++)
    {
        //cerr<<"Matching "<<i<<"\n";
        // single cluster stuff
        matchCluster(i);
        m_histos.FillHisto(50+hoffset, m_matchedClusters.size(), weight, sysNum);
        double responseMax = (m_matchedClusters[0]->calibratedEt() - m_genParticle->Et())/m_genParticle->Et();
        m_histos.FillHisto(55+hoffset, responseMax, weight, sysNum);
        for(const auto cluster : m_matchedClusters)
        {
            double deta = (cluster->eta()-m_genParticle->Eta());
            double dphi = TVector2::Phi_mpi_pi(cluster->phi()-m_genParticle->Phi());
            double response = (cluster->calibratedEt() - m_genParticle->Et())/m_genParticle->Et();
            m_histos.FillHisto(51+hoffset, deta, weight, sysNum);
            m_histos.FillHisto(52+hoffset, dphi, weight, sysNum);
            m_histos.FillHisto(53+hoffset, deta, dphi, weight, sysNum);
            m_histos.FillHisto(54+hoffset, response, weight, sysNum);

            // position with respect to main cluster
            if(cluster!=m_matchedClusters[0])
            {
                double detaSub = (cluster->eta()-m_matchedClusters[0]->eta());
                double dphiSub = TVector2::Phi_mpi_pi(cluster->phi()-m_matchedClusters[0]->phi());
                m_histos.FillHisto(61+hoffset, detaSub, weight, sysNum);
                m_histos.FillHisto(62+hoffset, dphiSub, weight, sysNum);
                m_histos.FillHisto(63+hoffset, detaSub, dphiSub, weight, sysNum);
            }
        }

        // supercluster stuff
        matchSuperCluster(i);
        if(!m_matchedSuperCluster) continue;
        //
        double deta = (m_matchedSuperCluster->eta()-m_genParticle->Eta());
        double dphi = TVector2::Phi_mpi_pi(m_matchedSuperCluster->phi()-m_genParticle->Phi());
        double response = (m_matchedSuperCluster->et() - m_genParticle->Et())/m_genParticle->Et();
        double reducedPhi = (m_genParticle->Phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        double localPhi = (reducedPhi - floor(reducedPhi))*(2*TMath::Pi()/18.);
        m_histos.FillHisto(70+hoffset, m_matchedSuperCluster->nClusters(), weight, sysNum);
        m_histos.FillHisto(71+hoffset, deta, weight, sysNum);
        m_histos.FillHisto(72+hoffset, dphi, weight, sysNum);
        m_histos.FillHisto(73+hoffset, response, weight, sysNum);
        m_histos.FillHisto(74+hoffset, response, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(78+hoffset, response, localPhi, weight, sysNum);
        if(response<-0.3)
        {
            m_histos.FillHisto(75+hoffset, fabs(m_genParticle->Eta()), m_genParticle->Phi(), weight, sysNum);
            m_histos.FillHisto(76+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            //m_histos.FillHisto(77+hoffset, m_genParticle->Phi(), weight, sysNum);
            m_histos.FillHisto(77+hoffset, localPhi, weight, sysNum);
        }


        // 
        matchSuperClusterHard(i);
        matchSuperClusterHardNoPUThreshold(i);
        const Tower* seed = m_matchedSuperCluster->cluster(0);
        const Tower* seedHard = (m_matchedSuperClusterHard ? findSeedClusterHard(m_matchedSuperCluster) : 0);
        const Tower* seedHardNoTh = (m_matchedSuperClusterHardNoPUThreshold ? findSeedClusterHardNoPUThreshold(m_matchedSuperCluster) : 0);
        //if(seedHard) cout<<"seed Et="<<m_matchedSuperCluster->cluster(0)->calibratedEt()<<" -> hard="<<seedHard->calibratedEt()<<"\n";
        //if(seedHardNoTh) cout<<"seed Et="<<m_matchedSuperCluster->cluster(0)->calibratedEt()<<" -> hard no th="<<seedHardNoTh->calibratedEt()<<"\n";
        int triggerRegion = m_egammaAlgo.triggerRegionIndex(seed->hits()[0]->eta(), seed->hits()[0]->zside(), seed->hits()[0]->sector(), seed->hits()[0]->subsector());
        int nhits = 0;
        nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 1);
        nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 2);
        nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 3);
        nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 4);
        nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 5);
        if(seedHard && seedHardNoTh)
        {
            m_histos.FillHisto(100+hoffset, seedHard->eta() - seed->eta(), weight, sysNum);
            m_histos.FillHisto(101+hoffset, TVector2::Phi_mpi_pi(seedHard->phi() - seed->phi()), weight, sysNum);
            m_histos.FillHisto(102+hoffset, seedHard->calibratedEnergy() - seed->calibratedEnergy(), weight, sysNum);
            m_histos.FillHisto(103+hoffset, seed->eta(), weight, sysNum);
            m_histos.FillHisto(104+hoffset, seed->calibratedEnergy(), weight, sysNum);
            m_histos.FillHisto(105+hoffset, nhits, weight, sysNum);
            m_histos.FillHisto(106+hoffset, seedHardNoTh->calibratedEnergy()/seedHard->calibratedEnergy(), weight, sysNum);

            m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);
        }
    }

}



/*****************************************************************/
void AnalysisClustering::matchCluster(int i)
/*****************************************************************/
{
    m_matchedClusters.clear();

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    for(const auto& cluster : m_clusters )
    {
        double deta = (cluster.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cluster.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.3)
        {
            m_matchedClusters.push_back(&cluster);
        }
    }
    sort(m_matchedClusters.begin(), m_matchedClusters.end(), AnHiMa::towerSort);

}

/*****************************************************************/
void AnalysisClustering::matchClusterHard(int i)
/*****************************************************************/
{
    m_matchedClustersHard.clear();

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    for(const auto& cluster : m_clustersHard )
    {
        double deta = (cluster.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cluster.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.3)
        {
            m_matchedClustersHard.push_back(&cluster);
        }
    }
    sort(m_matchedClustersHard.begin(), m_matchedClustersHard.end(), AnHiMa::towerSort);

}

/*****************************************************************/
void AnalysisClustering::matchClusterHardNoPUThreshold(int i)
/*****************************************************************/
{
    m_matchedClustersHardNoPUThreshold.clear();

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    for(const auto& cluster : m_clustersHardNoPUThreshold )
    {
        double deta = (cluster.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cluster.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.3)
        {
            m_matchedClustersHardNoPUThreshold.push_back(&cluster);
        }
    }
    sort(m_matchedClustersHardNoPUThreshold.begin(), m_matchedClustersHardNoPUThreshold.end(), AnHiMa::towerSort);

}

/*****************************************************************/
void AnalysisClustering::matchSuperCluster(int i)
/*****************************************************************/
{
    m_matchedSuperCluster = 0;

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClusters )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperCluster = &sc;
            maxEt = sc.et();
        }
    }
}

/*****************************************************************/
void AnalysisClustering::matchSuperClusterHard(int i)
/*****************************************************************/
{
    m_matchedSuperClusterHard = 0;

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClustersHard )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperClusterHard = &sc;
            maxEt = sc.et();
        }
    }
}

/*****************************************************************/
void AnalysisClustering::matchSuperClusterHardNoPUThreshold(int i)
/*****************************************************************/
{
    m_matchedSuperClusterHardNoPUThreshold = 0;

    m_genParticle = &event().genparticles()[i];
    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClustersHardNoPUThreshold )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperClusterHardNoPUThreshold = &sc;
            maxEt = sc.et();
        }
    }
}



/*****************************************************************/
const Tower* AnalysisClustering::findSeedClusterHard(const SuperCluster* sc)
/*****************************************************************/
{
    const Tower* seed = sc->cluster(0);
    const Tower* matchedSeed = 0;
    double eta = seed->eta();
    double phi = seed->phi();
    double minDR = 9999.;
    for(const auto& cl : m_matchedSuperClusterHard->clusters())
    {
        double deta = (cl->eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl->phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.05 && dr<minDR)
        {
            matchedSeed = cl;
            minDR = dr;
        }
    }
    return matchedSeed;
}

/*****************************************************************/
const Tower* AnalysisClustering::findSeedClusterHardNoPUThreshold(const SuperCluster* sc)
/*****************************************************************/
{
    const Tower* seed = sc->cluster(0);
    const Tower* matchedSeed = 0;
    double eta = seed->eta();
    double phi = seed->phi();
    double minDR = 9999.;
    for(const auto& cl : m_matchedSuperClusterHardNoPUThreshold->clusters())
    {
        double deta = (cl->eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl->phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.05 && dr<minDR)
        {
            matchedSeed = cl;
            minDR = dr;
        }
    }
    return matchedSeed;
}
