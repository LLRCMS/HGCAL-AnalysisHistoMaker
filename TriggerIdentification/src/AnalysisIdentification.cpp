/**
 *  @file  AnalysisIdentification.cpp
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

#include "AnHiMaHGCAL/TriggerIdentification/interface/AnalysisIdentification.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisIdentification::AnalysisIdentification():IAnalysis(),
    m_genParticle(0),
    m_superCluster(0),
    m_sample("signal"),
    m_etCut(10.)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisIdentification::~AnalysisIdentification()
/*****************************************************************/
{

}



/*****************************************************************/
bool AnalysisIdentification::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // Read parameters
    m_sample = m_reader.params().GetValue("Sample", "signal");
    m_etCut  = m_reader.params().GetValue("EtCut", 10.);

    m_egammaAlgo.initialize(event(), m_reader.params());



    return true;
}


/*****************************************************************/
void AnalysisIdentification::execute()
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

    if(m_sample=="signal")
    {
        for(int i=0;i<2;i++)
        {
            // get gen electron
            m_genParticle = &event().genparticles()[i];
            if(fabs(m_genParticle->Eta())<2.9 && fabs(m_genParticle->Eta())>1.6) 
            {
                // match super-cluster to gen electron
                matchSuperCluster();
                if(m_superCluster && m_superCluster->et()>m_etCut)
                {
                    fillHistos();
                }
            }
        }
    }
    else
    {
        for(const auto& sc : m_superClusters)
        {
            if(sc.et()>m_etCut)
            {
                m_superCluster = &sc;
                fillHistos();
            }
        }
    }
}

/*****************************************************************/
void AnalysisIdentification::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;



    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //
    double reducedPhiSC = (m_superCluster->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
    double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
    const Tower* seed = m_superCluster->cluster(0);
    m_histos.FillHisto(10+hoffset, seed->eta(), weight, sysNum); 
    m_histos.FillHisto(11+hoffset, seed->phi(), weight, sysNum); 
    m_histos.FillHisto(12+hoffset, seed->calibratedEnergy(), weight, sysNum); 
    m_histos.FillHisto(13+hoffset, seed->calibratedEt(), weight, sysNum); 
    m_histos.FillHisto(14+hoffset, localPhiSC, weight, sysNum); 
    //
    m_histos.FillHisto(20+hoffset, seed->layerCalibratedEnergy(20), weight, sysNum); 
    m_histos.FillHisto(21+hoffset, seed->layerCalibratedEnergy(21), weight, sysNum); 
    m_histos.FillHisto(22+hoffset, seed->layerCalibratedEnergy(22), weight, sysNum); 
    m_histos.FillHisto(23+hoffset, seed->layerCalibratedEnergy(23), weight, sysNum); 
    m_histos.FillHisto(24+hoffset, seed->layerCalibratedEnergy(24), weight, sysNum); 
    m_histos.FillHisto(25+hoffset, seed->layerCalibratedEnergy(25), weight, sysNum); 
    m_histos.FillHisto(26+hoffset, seed->layerCalibratedEnergy(26), weight, sysNum); 
    m_histos.FillHisto(27+hoffset, seed->layerCalibratedEnergy(27), weight, sysNum); 
    m_histos.FillHisto(28+hoffset, seed->layerCalibratedEnergy(28), weight, sysNum); 
    m_histos.FillHisto(29+hoffset, seed->layerCalibratedEnergy(29), weight, sysNum); 
    m_histos.FillHisto(30+hoffset, seed->layerCalibratedEnergy(30), weight, sysNum); 
    //
    m_histos.FillHisto(50+hoffset, seed->layerNHits(10), weight, sysNum); 
    m_histos.FillHisto(51+hoffset, seed->layerNHits(11), weight, sysNum); 
    m_histos.FillHisto(52+hoffset, seed->layerNHits(12), weight, sysNum); 
    m_histos.FillHisto(53+hoffset, seed->layerNHits(13), weight, sysNum); 
    m_histos.FillHisto(54+hoffset, seed->layerNHits(14), weight, sysNum); 
    m_histos.FillHisto(55+hoffset, seed->layerNHits(15), weight, sysNum); 
    m_histos.FillHisto(56+hoffset, seed->layerNHits(16), weight, sysNum); 
    m_histos.FillHisto(57+hoffset, seed->layerNHits(17), weight, sysNum); 
    m_histos.FillHisto(58+hoffset, seed->layerNHits(18), weight, sysNum); 
    m_histos.FillHisto(59+hoffset, seed->layerNHits(19), weight, sysNum); 
    m_histos.FillHisto(60+hoffset, seed->layerNHits(20), weight, sysNum);
    //
    int maxLayer = 0;
    int maxLayerH = 0;
    double maxLayerEnergy = 0.;
    int maxLayerHit = 0;
    int firstLayer = 0;
    int lastLayer = 0;
    for(unsigned l=1;l<=30;l++)
    {
        double energy = seed->layerCalibratedEnergy(l);
        int hits = seed->layerNHits(l);
        if(energy>maxLayerEnergy)
        {
            maxLayer = l;
            maxLayerEnergy = energy;
        }
        if(hits>maxLayerHit)
        {
            maxLayerH = l;
            maxLayerHit = hits;
        }
        if(firstLayer==0 && energy>0.) firstLayer = l;
    }
    for(unsigned l=30;l>=1;l--)
    {
        double energy = seed->layerCalibratedEnergy(l);
        if(lastLayer==0 && energy>0.) 
        {
            lastLayer = l;
            break;
        }
    }
    m_histos.FillHisto(100+hoffset, maxLayer, weight, sysNum);
    m_histos.FillHisto(101+hoffset, firstLayer, weight, sysNum);
    m_histos.FillHisto(102+hoffset, lastLayer, weight, sysNum);
    m_histos.FillHisto(103+hoffset, maxLayerH, weight, sysNum);

    /// BDT
    double bdt = m_egammaAlgo.bdtOutput(*seed);
    m_histos.FillHisto(120+hoffset, bdt, weight, sysNum);

    m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);
}





/*****************************************************************/
void AnalysisIdentification::matchSuperCluster()
/*****************************************************************/
{
    m_superCluster = 0;

    double eta = m_genParticle->Eta();
    double phi = m_genParticle->Phi();
    double maxEt = 0.;
    for(const auto& sc : m_superClusters )
    {
        double deta = (sc.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(sc.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && sc.et()>maxEt && (sc.et()-m_genParticle->Pt())/m_genParticle->Pt()>-0.5)
        {
            m_superCluster = &sc;
            maxEt = sc.et();
        }
    }
}







