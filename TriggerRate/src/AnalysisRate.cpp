/**
 *  @file  AnalysisRate.cpp
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

#include "AnHiMaHGCAL/TriggerRate/interface/AnalysisRate.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisRate::AnalysisRate():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisRate::~AnalysisRate()
/*****************************************************************/
{

}



/*****************************************************************/
bool AnalysisRate::initialize(const string& parameterFile)
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
void AnalysisRate::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;

    m_seeds.clear();
    m_clusters.clear();
    m_superClusters.clear();
    m_sortedSuperClusters.clear();
    m_sortedSuperClustersLooseID.clear();
    m_sortedSuperClustersMediumID.clear();

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
    // sort superclusters
    for(const auto& sc : m_superClusters) m_sortedSuperClusters.push_back(&sc);
    sort(m_sortedSuperClusters.begin(), m_sortedSuperClusters.end(), AnHiMa::superClusterSort);
    // filter identified clusters
    applyIdentification();

    fillHistos();
}

/*****************************************************************/
void AnalysisRate::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    const SuperCluster* maxCandidate = (m_sortedSuperClusters.size()>0 ? m_sortedSuperClusters[0] : 0);
    const SuperCluster* maxCandidateLooseID = (m_sortedSuperClustersLooseID.size()>0 ? m_sortedSuperClustersLooseID[0] : 0);
    const SuperCluster* maxCandidateMediumID = (m_sortedSuperClustersMediumID.size()>0 ? m_sortedSuperClustersMediumID[0] : 0);

    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //

    m_histos.FillHisto(100+hoffset, (maxCandidate ? maxCandidate->et() : 0.) , weight, sysNum);
    m_histos.FillHisto(101+hoffset, (maxCandidate ? fabs(maxCandidate->eta()) : 0.) , weight, sysNum);
    //
    m_histos.FillHisto(110+hoffset, (maxCandidateLooseID ? maxCandidateLooseID->et() : 0.) , weight, sysNum);
    m_histos.FillHisto(111+hoffset, (maxCandidateLooseID ? fabs(maxCandidateLooseID->eta()) : 0.) , weight, sysNum);
    //
    m_histos.FillHisto(120+hoffset, (maxCandidateMediumID ? maxCandidateMediumID->et() : 0.) , weight, sysNum);
    m_histos.FillHisto(121+hoffset, (maxCandidateMediumID ? fabs(maxCandidateMediumID->eta()) : 0.) , weight, sysNum);
}



/*****************************************************************/
void AnalysisRate::applyIdentification()
/*****************************************************************/
{
    for(const auto sc : m_sortedSuperClusters)
    {
        double bdt = m_egammaAlgo.bdtOutput(*sc->cluster(0));
        bool passLoose = (bdt>=-0.84);
        bool passMedium = (bdt>=-0.76);
        double reducedPhiSC = (sc->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
        bool regionLoose = (fabs(sc->eta())<2.2 || sc->et()<30. || localPhiSC<0.05 || localPhiSC>(2*TMath::Pi()/18.-0.05));
        if((regionLoose && passLoose) || (!regionLoose && passMedium))
        {
            m_sortedSuperClustersLooseID.push_back(sc);
        }
        if(passMedium)
        {
            m_sortedSuperClustersMediumID.push_back(sc);
        }
    }
}




