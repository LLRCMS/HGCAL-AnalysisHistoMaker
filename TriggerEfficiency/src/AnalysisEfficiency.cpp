/**
 *  @file  AnalysisEfficiency.cpp
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

#include "AnHiMaHGCAL/TriggerEfficiency/interface/AnalysisEfficiency.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisEfficiency::AnalysisEfficiency():IAnalysis(),
    m_genParticle(0),
    m_matchedSeed(0),
    m_matchedCluster(0),
    m_matchedSuperCluster(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisEfficiency::~AnalysisEfficiency()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisEfficiency::initialize(const string& parameterFile)
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
void AnalysisEfficiency::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_clusters.clear();
    m_superClusters.clear();

    // seeding
    m_egammaAlgo.seeding(event(), m_seeds);
    // clustering
    m_egammaAlgo.clustering(event(), m_seeds, m_clusters);
    // cluster energy corrections
    m_egammaAlgo.clusterPileupSubtraction(m_clusters);
    m_egammaAlgo.clusterThresholdCorrection(m_clusters);
    // super-clustering
    m_egammaAlgo.superClustering(m_clusters, m_superClusters);
    // supercluster energy corrections
    m_egammaAlgo.superClusterCorrection(m_superClusters);

    for(int i=0;i<2;i++)
    {
        // get gen electron
        m_genParticle = &event().genparticles()[i];
        if(fabs(m_genParticle->Eta())<2.9 && fabs(m_genParticle->Eta())>1.6) 
        {
            // match seed, cluster and super-cluster to gen electron
            matchSeed();
            matchCluster();
            matchSuperCluster();
            fillHistos();
        }
    }
}

/*****************************************************************/
void AnalysisEfficiency::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;



    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //
    double reducedPhi = (m_genParticle->Phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
    double localPhi = (reducedPhi - floor(reducedPhi))*(2*TMath::Pi()/18.);
    if(m_genParticle)
    {

        m_histos.FillHisto(10+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(11+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(12+hoffset, event().npu(), weight, sysNum);
    }
    //
    if(m_matchedSeed)
    {
        m_histos.FillHisto(20+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(21+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(22+hoffset, event().npu(), weight, sysNum);
    }
    //
    if(m_matchedCluster && m_matchedCluster->calibratedEt()>5.)
    {
        m_histos.FillHisto(30+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(31+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(32+hoffset, event().npu(), weight, sysNum);
    }
    //
    if(m_matchedSuperCluster && m_matchedSuperCluster->et()>5.)
    {
        m_histos.FillHisto(40+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(41+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(42+hoffset, event().npu(), weight, sysNum);
    }

}


/*****************************************************************/
void AnalysisEfficiency::matchSeed()
/*****************************************************************/
{
    m_matchedSeed = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& seed : m_seeds)
    {
        double deta = (seed.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(seed.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && seed.calibratedEt()>maxEt)
        {
            m_matchedSeed = &seed;
            maxEt = seed.calibratedEt();
        }
    }
}



/*****************************************************************/
void AnalysisEfficiency::matchCluster()
/*****************************************************************/
{
    m_matchedCluster = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& cluster : m_clusters )
    {
        double deta = (cluster.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cluster.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && cluster.calibratedEt()>maxEt)
        {
            m_matchedCluster = &cluster;
            maxEt = cluster.calibratedEt();
        }
    }
}



/*****************************************************************/
void AnalysisEfficiency::matchSuperCluster()
/*****************************************************************/
{
    m_matchedSuperCluster = 0;

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







