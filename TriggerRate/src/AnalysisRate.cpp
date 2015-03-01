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


    /// get BDT cuts
    string fileNameBDTcuts = m_reader.params().GetValue("FileBDTCuts", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/bdtCutsVsEta_QCD.root");
    TFile* fileBDTcuts = TFile::Open(fileNameBDTcuts.c_str());
    TGraph* bdtCuts995 = (TGraph*)fileBDTcuts->Get("bdtCutVsEta_eff0.995");
    TGraph* bdtCuts99  = (TGraph*)fileBDTcuts->Get("bdtCutVsEta_eff0.99");
    TGraph* bdtCuts985 = (TGraph*)fileBDTcuts->Get("bdtCutVsEta_eff0.985");
    TGraph* bdtCuts98  = (TGraph*)fileBDTcuts->Get("bdtCutVsEta_eff0.98");
    TGraph* bdtCuts975 = (TGraph*)fileBDTcuts->Get("bdtCutVsEta_eff0.975");
    TGraph* bdtCuts97  = (TGraph*)fileBDTcuts->Get("bdtCutVsEta_eff0.97");
    m_bdtCuts[995] = bdtCuts995;
    m_bdtCuts[99]  = bdtCuts99;
    m_bdtCuts[985] = bdtCuts985;
    m_bdtCuts[98]  = bdtCuts98;
    m_bdtCuts[975] = bdtCuts975;
    m_bdtCuts[97]  = bdtCuts97;
    for(auto& eff_graph : m_bdtCuts)
    {
        TGraph* g = eff_graph.second;
        assert(g);
    }
    fileBDTcuts->Close();


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
    //
    const SuperCluster* max2ndCandidate = (m_sortedSuperClusters.size()>1 ? m_sortedSuperClusters[1] : 0);
    const SuperCluster* max2ndCandidateLooseID = (m_sortedSuperClustersLooseID.size()>1 ? m_sortedSuperClustersLooseID[1] : 0);
    const SuperCluster* max2ndCandidateMediumID = (m_sortedSuperClustersMediumID.size()>1 ? m_sortedSuperClustersMediumID[1] : 0);

    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    m_histos.FillHisto(10+hoffset, m_seeds.size(), weight, sysNum);
    //

    //m_histos.FillHisto(100+hoffset, (maxCandidate && maxCandidate->et()<100. ? maxCandidate->et() : 0.) , weight, sysNum);
    //m_histos.FillHisto(101+hoffset, (maxCandidate && maxCandidate->et()<100. ? fabs(maxCandidate->eta()) : 0.) , weight, sysNum);
    //m_histos.FillHisto(102+hoffset, (max2ndCandidate && max2ndCandidate->et()<100. ? max2ndCandidate->et() : 0.) , weight, sysNum);
    ////
    //m_histos.FillHisto(110+hoffset, (maxCandidateLooseID && maxCandidateLooseID->et()<100. ? maxCandidateLooseID->et() : 0.) , weight, sysNum);
    //m_histos.FillHisto(111+hoffset, (maxCandidateLooseID && maxCandidateLooseID->et()<100. ? fabs(maxCandidateLooseID->eta()) : 0.) , weight, sysNum);
    //m_histos.FillHisto(112+hoffset, (max2ndCandidateLooseID && max2ndCandidateLooseID->et()<100. ? max2ndCandidateLooseID->et() : 0.) , weight, sysNum);
    ////
    //m_histos.FillHisto(120+hoffset, (maxCandidateMediumID && maxCandidateMediumID->et()<100. ? maxCandidateMediumID->et() : 0.) , weight, sysNum);
    //m_histos.FillHisto(121+hoffset, (maxCandidateMediumID && maxCandidateMediumID->et()<100. ? fabs(maxCandidateMediumID->eta()) : 0.) , weight, sysNum);
    //m_histos.FillHisto(122+hoffset, (max2ndCandidateMediumID && max2ndCandidateMediumID->et()<100. ? max2ndCandidateMediumID->et() : 0.) , weight, sysNum);

    m_histos.FillHisto(100+hoffset, (maxCandidate  ? maxCandidate->et() : 0.) , weight, sysNum);
    m_histos.FillHisto(101+hoffset, (maxCandidate  ? fabs(maxCandidate->eta()) : 0.) , weight, sysNum);
    m_histos.FillHisto(102+hoffset, (max2ndCandidate  ? max2ndCandidate->et() : 0.) , weight, sysNum);
    //
    m_histos.FillHisto(110+hoffset, (maxCandidateLooseID  ? maxCandidateLooseID->et() : 0.) , weight, sysNum);
    m_histos.FillHisto(111+hoffset, (maxCandidateLooseID  ? fabs(maxCandidateLooseID->eta()) : 0.) , weight, sysNum);
    m_histos.FillHisto(112+hoffset, (max2ndCandidateLooseID  ? max2ndCandidateLooseID->et() : 0.) , weight, sysNum);
    //
    m_histos.FillHisto(120+hoffset, (maxCandidateMediumID  ? maxCandidateMediumID->et() : 0.) , weight, sysNum);
    m_histos.FillHisto(121+hoffset, (maxCandidateMediumID  ? fabs(maxCandidateMediumID->eta()) : 0.) , weight, sysNum);
    m_histos.FillHisto(122+hoffset, (max2ndCandidateMediumID  ? max2ndCandidateMediumID->et() : 0.) , weight, sysNum);

    for(const auto sc : m_sortedSuperClusters)
    {
        if(sc->et()<5. /*|| sc->et()>100.*/) continue;
        double reducedPhiSC = (sc->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
        m_histos.FillHisto(200+hoffset, sc->et() , weight, sysNum);
        m_histos.FillHisto(201+hoffset, fabs(sc->eta()) , weight, sysNum);
        m_histos.FillHisto(202+hoffset, localPhiSC , weight, sysNum);
    }
    for(const auto sc : m_sortedSuperClustersLooseID)
    {
        if(sc->et()<5. /*|| sc->et()>100.*/) continue;
        double reducedPhiSC = (sc->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
        m_histos.FillHisto(210+hoffset, sc->et() , weight, sysNum);
        m_histos.FillHisto(211+hoffset, fabs(sc->eta()) , weight, sysNum);
        m_histos.FillHisto(212+hoffset, localPhiSC , weight, sysNum);
    }
    for(const auto sc : m_sortedSuperClustersMediumID)
    {
        if(sc->et()<5. /*|| sc->et()>100.*/) continue;
        double reducedPhiSC = (sc->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
        m_histos.FillHisto(220+hoffset, sc->et() , weight, sysNum);
        m_histos.FillHisto(221+hoffset, fabs(sc->eta()) , weight, sysNum);
        m_histos.FillHisto(222+hoffset, localPhiSC , weight, sysNum);
    }
}



/*****************************************************************/
void AnalysisRate::applyIdentification()
/*****************************************************************/
{
    for(const auto sc : m_sortedSuperClusters)
    {
        // FIXME: create several lists with several BDTs
        string bdtName = "BDTG";
        if(m_egammaAlgo.useHalfLayers()) bdtName = "BDTG_MinBias_HalfLayers";
        double bdt = m_egammaAlgo.bdtOutput(*sc->cluster(0), bdtName.c_str());
        //bool pass995 = (bdt>=m_bdtCuts[995]->Eval(fabs(sc->eta())));
        bool pass99  = (bdt>=m_bdtCuts[99] ->Eval(fabs(sc->eta())));
        bool pass985 = (bdt>=m_bdtCuts[985]->Eval(fabs(sc->eta())));
        bool pass98  = (bdt>=m_bdtCuts[98] ->Eval(fabs(sc->eta())));
        //bool pass975 = (bdt>=m_bdtCuts[975]->Eval(fabs(sc->eta())));
        //bool pass97  = (bdt>=m_bdtCuts[97] ->Eval(fabs(sc->eta())));
        //double reducedPhiSC = (sc->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        //double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
        //bool regionLoose = (fabs(sc->eta())<2.2 || sc->et()<30. || localPhiSC<0.05 || localPhiSC>(2*TMath::Pi()/18.-0.05));
        if((sc->et()<30. && pass99) || (sc->et()>30. && pass98))
        {
            m_sortedSuperClustersLooseID.push_back(sc);
        }
        if(pass985)
        {
            m_sortedSuperClustersMediumID.push_back(sc);
        }
    }
}




