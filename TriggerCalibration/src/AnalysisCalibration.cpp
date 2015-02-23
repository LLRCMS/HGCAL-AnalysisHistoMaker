/**
 *  @file  AnalysisCalibration.cpp
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

#include "AnHiMaHGCAL/TriggerCalibration/interface/AnalysisCalibration.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisCalibration::AnalysisCalibration():IAnalysis(),
    m_genParticle(0),
    m_matchedSuperCluster(0),
    m_matchedSuperClusterHardNoPUThreshold(0),
    m_matchedSuperClusterCorrPU(0),
    m_matchedSuperClusterCorrTh(0),
    m_matchedSuperClusterCorrFull(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisCalibration::~AnalysisCalibration()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisCalibration::initialize(const string& parameterFile)
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
void AnalysisCalibration::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_clusters.clear();
    m_clustersHard.clear();
    m_clustersHardNoPUThreshold.clear();
    m_clustersCorrPU.clear();
    m_clustersCorrTh.clear();
    m_superClusters.clear();
    m_superClustersHard.clear();
    m_superClustersHardNoPUThreshold.clear();
    m_superClustersCorrPU.clear();
    m_superClustersCorrTh.clear();
    m_superClustersCorrFull.clear();
    // seeding
    m_egammaAlgo.seeding(event(), m_seeds);
    // clustering
    m_egammaAlgo.clustering(event(), m_seeds, m_clusters);
    m_egammaAlgo.clustering(event(), m_seeds, m_clustersHard, true, true);
    m_egammaAlgo.clustering(event(), m_seeds, m_clustersHardNoPUThreshold, true, false);
    // cluster energy corrections
    m_clustersCorrPU = m_clusters;
    m_egammaAlgo.clusterPileupSubtraction(m_clustersCorrPU);
    m_clustersCorrTh = m_clustersCorrPU;
    m_egammaAlgo.clusterThresholdCorrection(m_clustersCorrTh);
    // superClustering
    m_egammaAlgo.superClustering(m_clusters, m_superClusters);
    m_egammaAlgo.superClustering(m_clustersHard, m_superClustersHard);
    m_egammaAlgo.superClustering(m_clustersHardNoPUThreshold, m_superClustersHardNoPUThreshold);
    m_egammaAlgo.superClustering(m_clustersCorrPU, m_superClustersCorrPU);
    m_egammaAlgo.superClustering(m_clustersCorrTh, m_superClustersCorrTh);
    // supercluster energy corrections
    m_superClustersCorrFull = m_superClustersCorrTh;
    m_egammaAlgo.superClusterCorrection(m_superClustersCorrFull);


    for(int i=0;i<2;i++)
    {
        // get gen electron
        m_genParticle = &event().genparticles()[i];
        if(fabs(m_genParticle->Eta())<2.9 && fabs(m_genParticle->Eta())>1.6) 
        {
            matchSuperCluster();
            matchSuperClusterHard();
            matchSuperClusterHardNoPUThreshold();
            matchSuperClusterCorrPU();
            matchSuperClusterCorrTh();
            matchSuperClusterCorrFull();
            // ideal clustering
            m_idealClusters.clear();
            m_egammaAlgo.idealClustering(event(), m_idealClusters, m_genParticle->Eta(), m_genParticle->Phi());
            m_matchedIdealCluster = (m_idealClusters.size()>0 ? &m_idealClusters[0] : 0);
            if(m_matchedSuperCluster && m_matchedSuperClusterHard && m_matchedSuperClusterHardNoPUThreshold && m_matchedIdealCluster)
            {
                fillHistos();
            }
        }
    }
}

/*****************************************************************/
void AnalysisCalibration::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;



    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    // compare with generated
    double responseRaw      = (m_matchedSuperCluster->et() - m_genParticle->Pt())/m_genParticle->Pt();
    double responseHard     = (m_matchedSuperClusterHard->et() - m_genParticle->Pt())/m_genParticle->Pt();
    double responseHardNoTh = (m_matchedSuperClusterHardNoPUThreshold->et() - m_genParticle->Pt())/m_genParticle->Pt();
    double responseCorrPU   = (m_matchedSuperClusterCorrPU->et() - m_genParticle->Pt())/m_genParticle->Pt();
    double responseCorrTh   = (m_matchedSuperClusterCorrTh->et() - m_genParticle->Pt())/m_genParticle->Pt();
    double responseCorrFull = (m_matchedSuperClusterCorrFull->et() - m_genParticle->Pt())/m_genParticle->Pt();
    double responseIdeal    = (m_matchedIdealCluster->calibratedEt() - m_genParticle->Pt())/m_genParticle->Pt();

    m_histos.FillHisto(10+hoffset, responseRaw     , weight, sysNum);
    m_histos.FillHisto(11+hoffset, responseHard    , weight, sysNum);
    m_histos.FillHisto(12+hoffset, responseHardNoTh, weight, sysNum);
    m_histos.FillHisto(13+hoffset, responseCorrPU  , weight, sysNum);
    m_histos.FillHisto(14+hoffset, responseCorrTh  , weight, sysNum);
    m_histos.FillHisto(15+hoffset, responseCorrFull, weight, sysNum);
    m_histos.FillHisto(16+hoffset, responseIdeal   , weight, sysNum);

    // impact of individual corrections
    const Tower* seedRaw      = m_matchedSuperCluster->cluster(0);
    const Tower* seedHard     = findSeedClusterHard(m_matchedSuperCluster);
    const Tower* seedHardNoTh = findSeedClusterHardNoPUThreshold(m_matchedSuperCluster);
    const Tower* seedCorrPU   = m_matchedSuperClusterCorrPU->cluster(0);
    const Tower* seedCorrTh   = m_matchedSuperClusterCorrTh->cluster(0);

    double diffPUraw        = (seedRaw->calibratedEt() - seedHard->calibratedEt());
    double diffPUcorr       = (seedCorrPU->calibratedEt() - seedHard->calibratedEt());
    double responseThNoPU   = (seedHard->calibratedEt() - seedHardNoTh->calibratedEt())/seedHardNoTh->calibratedEt();
    double responseThCorrPU = (seedCorrPU->calibratedEt() - seedHardNoTh->calibratedEt())/seedHardNoTh->calibratedEt();
    double responseThCorrTh = (seedCorrTh->calibratedEt() - seedHardNoTh->calibratedEt())/seedHardNoTh->calibratedEt();

    m_histos.FillHisto(50+hoffset, diffPUraw       , weight, sysNum);
    m_histos.FillHisto(51+hoffset, diffPUcorr      , weight, sysNum);
    m_histos.FillHisto(52+hoffset, responseThNoPU  , weight, sysNum);
    m_histos.FillHisto(53+hoffset, responseThCorrPU, weight, sysNum);
    m_histos.FillHisto(54+hoffset, responseThCorrTh, weight, sysNum);


    // layer energy fractions in ideal clusters
    for(int l=1;l<=30;l++)   
    {
        double layerEnergy = m_matchedIdealCluster->layerEnergy(l);
        m_histos.Fill1BinHisto(60+hoffset, l, layerEnergy/m_matchedIdealCluster->energy(), weight, sysNum);
    }

    // cluster corrections
    double corrPU = (seedHard->calibratedEnergy() - seedRaw->calibratedEnergy())/(double)seedRaw->nHits();
    double corrTh = seedHardNoTh->calibratedEnergy()/seedHard->calibratedEnergy();

    int triggerRegion = m_egammaAlgo.triggerRegionIndex(seedRaw->hits()[0]->eta(), seedRaw->hits()[0]->zside(), seedRaw->hits()[0]->sector(), seedRaw->hits()[0]->subsector());
    int nhits = 0;
    bool useHalfLayers = m_reader.params().GetValue("UseHalfLayers", false);
    if(!useHalfLayers) nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 1);
    nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 2);
    if(!useHalfLayers) nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 3);
    nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 4);
    if(!useHalfLayers) nhits += m_egammaAlgo.triggerRegionHits(triggerRegion, 5);

    m_histos.FillHisto(100+hoffset, corrPU, weight, sysNum);
    m_histos.FillHisto(101+hoffset, corrTh, weight, sysNum);
    m_histos.FillHisto(102+hoffset, seedRaw->eta(), weight, sysNum);
    m_histos.FillHisto(103+hoffset, nhits, weight, sysNum);
    m_histos.FillHisto(104+hoffset, seedRaw->calibratedEnergy(), weight, sysNum);
    m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);

    // supercluster corrections
    if(m_matchedSuperClusterCorrTh->et()>7.)
    {
        double corrGen = m_genParticle->Pt()/m_matchedSuperClusterCorrTh->et();

        double scReducedPhi = (m_matchedSuperClusterCorrTh->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        double scLocalPhi = (scReducedPhi - floor(scReducedPhi))*(2*TMath::Pi()/18.);

        m_histos.FillHisto(110+hoffset, corrGen, weight, sysNum);
        m_histos.FillHisto(111+hoffset, m_matchedSuperClusterCorrTh->eta(), weight, sysNum);
        m_histos.FillHisto(112+hoffset, scLocalPhi, weight, sysNum);
        m_histos.FillHisto(113+hoffset, m_matchedSuperClusterCorrTh->energy(), weight, sysNum);
        m_histos.FillHisto(114+hoffset, m_matchedSuperClusterCorrTh->et(), weight, sysNum);
        m_histos.FillHisto(115+hoffset, m_matchedSuperClusterCorrTh->nClusters(), weight, sysNum);
        m_histos.FillNtuple(600+hoffset, event().run(), event().event(), weight, sysNum);
    }



}




/*****************************************************************/
void AnalysisCalibration::matchSuperCluster()
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

/*****************************************************************/
void AnalysisCalibration::matchSuperClusterHard()
/*****************************************************************/
{
    m_matchedSuperClusterHard = 0;

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
void AnalysisCalibration::matchSuperClusterHardNoPUThreshold()
/*****************************************************************/
{
    m_matchedSuperClusterHardNoPUThreshold = 0;

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
void AnalysisCalibration::matchSuperClusterCorrPU()
/*****************************************************************/
{
    m_matchedSuperClusterCorrPU = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClustersCorrPU )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperClusterCorrPU = &sc;
            maxEt = sc.et();
        }
    }
}

/*****************************************************************/
void AnalysisCalibration::matchSuperClusterCorrTh()
/*****************************************************************/
{
    m_matchedSuperClusterCorrTh = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClustersCorrTh )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperClusterCorrTh = &sc;
            maxEt = sc.et();
        }
    }
}

/*****************************************************************/
void AnalysisCalibration::matchSuperClusterCorrFull()
/*****************************************************************/
{
    m_matchedSuperClusterCorrFull = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClustersCorrFull )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt)
        {
            m_matchedSuperClusterCorrFull = &sc;
            maxEt = sc.et();
        }
    }
}


/*****************************************************************/
const Tower* AnalysisCalibration::findSeedClusterHard(const SuperCluster* sc)
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
const Tower* AnalysisCalibration::findSeedClusterHardNoPUThreshold(const SuperCluster* sc)
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
