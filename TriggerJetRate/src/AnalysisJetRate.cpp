/**
 *  @file  AnalysisJetRate.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    21/02/2015
 *
 *  @internal
 *     Created :  21/02/2015
 * Last update :  21/02/2015 13:45:14
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/TriggerJetRate/interface/AnalysisJetRate.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisJetRate::AnalysisJetRate():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisJetRate::~AnalysisJetRate()
/*****************************************************************/
{

}



/*****************************************************************/
bool AnalysisJetRate::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_egammaAlgo.initialize(event(), m_reader.params());
    m_jetAlgo.initialize(event(), m_reader.params());


    m_jets["Size01"] = vector<Tower>();
    m_jets["Size02"] = vector<Tower>();

    return true;
}


/*****************************************************************/
void AnalysisJetRate::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;


    m_emSeeds.clear();
    m_emClusters.clear();
    m_hadClusters.clear();
    m_clusters.clear();
    m_superClusters.clear();
    for(auto& col_jets : m_jets) col_jets.second.clear();
    m_sortedJets.clear();

    //-------------------------------------------
    // e/g algo
    m_egammaAlgo.seeding(event(), m_emSeeds);
    m_egammaAlgo.clustering(event(), m_emSeeds, m_emClusters);
    m_egammaAlgo.clusterPileupSubtraction(m_emClusters);
    m_egammaAlgo.clusterThresholdCorrection(m_emClusters);
    //-------------------------------------------
    // Jet algo
    m_jetAlgo.hcalSeeding(event(), event().simhitsFH(), m_hadClusters);
    //
    m_clusters = m_emClusters;
    for(const auto cl : m_hadClusters) if(cl.calibratedEt()>1. && cl.nHits()>5) m_clusters.push_back(cl);
    m_jetAlgo.superClustering(event(), m_clusters, m_superClusters, 0.4);
    //
    vector<SimHit> simhits = event().simhits();
    simhits.insert(simhits.end(), event().simhitsFH().begin(), event().simhitsFH().end());
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets["Size01"], 0.1); 
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets["Size02"], 0.2); 
    //
    m_jetAlgo.pileupSubtraction(m_jets["Size01"], "0.1");
    m_jetAlgo.pileupSubtraction(m_jets["Size02"], "0.2");
    m_jetAlgo.thresholdCorrection(m_jets["Size01"], "0.1");
    m_jetAlgo.thresholdCorrection(m_jets["Size02"], "0.2");
    m_jetAlgo.truthCorrection(m_jets["Size02"], m_jets["Size01"]);

    // sort jets
    for(const auto& jet : m_jets["Size02"])
    {
        if(jet.calibratedEt()>0.) m_sortedJets.push_back(&jet); 
    }
    sort(m_sortedJets.begin(), m_sortedJets.end(), AnHiMa::towerSort);

    fillHistos();
}

/*****************************************************************/
void AnalysisJetRate::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    const Tower* maxCandidate = (m_sortedJets.size()>0 ? m_sortedJets[0] : 0);
    //
    const Tower* max2ndCandidate = (m_sortedJets.size()>1 ? m_sortedJets[1] : 0);


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events
    //



    m_histos.FillHisto(100+hoffset, (maxCandidate  ? maxCandidate->calibratedEt() : 0.) , weight, sysNum);
    m_histos.FillHisto(101+hoffset, (maxCandidate  ? fabs(maxCandidate->eta()) : 0.) , weight, sysNum);
    m_histos.FillHisto(102+hoffset, (max2ndCandidate  ? max2ndCandidate->calibratedEt() : 0.) , weight, sysNum);
    //
    bool foundOppositeEtaJets = false;
    for(unsigned i=0;i<m_sortedJets.size();i++)
    {
        const auto& jet1 = m_sortedJets.at(i);
        for(unsigned j=i+1;j<m_sortedJets.size();j++)
        {
            const auto& jet2 = m_sortedJets.at(j);
            if(jet1->eta()*jet2->eta()<0.)
            {
                m_histos.FillHisto(103+hoffset, jet2->calibratedEt(), weight, sysNum);
                foundOppositeEtaJets = true;
                break;
            }
        }
        if(foundOppositeEtaJets) break;
    }
    if(!foundOppositeEtaJets) m_histos.FillHisto(103+hoffset, 0., weight, sysNum);
    // M > 200 GeV
    bool foundOppositeEtaJetsM200 = false;
    for(unsigned i=0;i<m_sortedJets.size();i++)
    {
        const auto& jet1 = m_sortedJets.at(i);
        TLorentzVector v1;
        v1.SetPtEtaPhiM(jet1->calibratedEt(), jet1->eta(), jet1->phi(), 0.);
        for(unsigned j=i+1;j<m_sortedJets.size();j++)
        {
            const auto& jet2 = m_sortedJets.at(j);
            TLorentzVector v2;
            v2.SetPtEtaPhiM(jet2->calibratedEt(), jet2->eta(), jet2->phi(), 0.);
            TLorentzVector v1v2(v2);
            v1v2 += v1;
            if(jet1->eta()*jet2->eta()<0. && v1v2.M()>200.)
            {
                m_histos.FillHisto(104+hoffset, jet2->calibratedEt(), weight, sysNum);
                foundOppositeEtaJetsM200 = true;
                break;
            }
        }
        if(foundOppositeEtaJetsM200) break;
    }
    if(!foundOppositeEtaJetsM200) m_histos.FillHisto(104+hoffset, 0., weight, sysNum);
    // M > 250 GeV
    bool foundOppositeEtaJetsM250 = false;
    for(unsigned i=0;i<m_sortedJets.size();i++)
    {
        const auto& jet1 = m_sortedJets.at(i);
        TLorentzVector v1;
        v1.SetPtEtaPhiM(jet1->calibratedEt(), jet1->eta(), jet1->phi(), 0.);
        for(unsigned j=i+1;j<m_sortedJets.size();j++)
        {
            const auto& jet2 = m_sortedJets.at(j);
            TLorentzVector v2;
            v2.SetPtEtaPhiM(jet2->calibratedEt(), jet2->eta(), jet2->phi(), 0.);
            TLorentzVector v1v2(v2);
            v1v2 += v1;
            if(jet1->eta()*jet2->eta()<0. && v1v2.M()>250.)
            {
                m_histos.FillHisto(105+hoffset, jet2->calibratedEt(), weight, sysNum);
                foundOppositeEtaJetsM250 = true;
                break;
            }
        }
        if(foundOppositeEtaJetsM250) break;
    }
    if(!foundOppositeEtaJetsM250) m_histos.FillHisto(105+hoffset, 0., weight, sysNum);
    // M > 300 GeV
    bool foundOppositeEtaJetsM300 = false;
    for(unsigned i=0;i<m_sortedJets.size();i++)
    {
        const auto& jet1 = m_sortedJets.at(i);
        TLorentzVector v1;
        v1.SetPtEtaPhiM(jet1->calibratedEt(), jet1->eta(), jet1->phi(), 0.);
        for(unsigned j=i+1;j<m_sortedJets.size();j++)
        {
            const auto& jet2 = m_sortedJets.at(j);
            TLorentzVector v2;
            v2.SetPtEtaPhiM(jet2->calibratedEt(), jet2->eta(), jet2->phi(), 0.);
            TLorentzVector v1v2(v2);
            v1v2 += v1;
            if(jet1->eta()*jet2->eta()<0. && v1v2.M()>300.)
            {
                m_histos.FillHisto(106+hoffset, jet2->calibratedEt(), weight, sysNum);
                foundOppositeEtaJetsM300 = true;
                break;
            }
        }
        if(foundOppositeEtaJetsM300) break;
    }
    if(!foundOppositeEtaJetsM300) m_histos.FillHisto(106+hoffset, 0., weight, sysNum);
    // M > 400 GeV
    bool foundOppositeEtaJetsM400 = false;
    for(unsigned i=0;i<m_sortedJets.size();i++)
    {
        const auto& jet1 = m_sortedJets.at(i);
        TLorentzVector v1;
        v1.SetPtEtaPhiM(jet1->calibratedEt(), jet1->eta(), jet1->phi(), 0.);
        for(unsigned j=i+1;j<m_sortedJets.size();j++)
        {
            const auto& jet2 = m_sortedJets.at(j);
            TLorentzVector v2;
            v2.SetPtEtaPhiM(jet2->calibratedEt(), jet2->eta(), jet2->phi(), 0.);
            TLorentzVector v1v2(v2);
            v1v2 += v1;
            if(jet1->eta()*jet2->eta()<0. && v1v2.M()>400.)
            {
                m_histos.FillHisto(107+hoffset, jet2->calibratedEt(), weight, sysNum);
                foundOppositeEtaJetsM400 = true;
                break;
            }
        }
        if(foundOppositeEtaJetsM400) break;
    }
    if(!foundOppositeEtaJetsM400) m_histos.FillHisto(107+hoffset, 0., weight, sysNum);
    // M > 500 GeV
    bool foundOppositeEtaJetsM500 = false;
    for(unsigned i=0;i<m_sortedJets.size();i++)
    {
        const auto& jet1 = m_sortedJets.at(i);
        TLorentzVector v1;
        v1.SetPtEtaPhiM(jet1->calibratedEt(), jet1->eta(), jet1->phi(), 0.);
        for(unsigned j=i+1;j<m_sortedJets.size();j++)
        {
            const auto& jet2 = m_sortedJets.at(j);
            TLorentzVector v2;
            v2.SetPtEtaPhiM(jet2->calibratedEt(), jet2->eta(), jet2->phi(), 0.);
            TLorentzVector v1v2(v2);
            v1v2 += v1;
            if(jet1->eta()*jet2->eta()<0. && v1v2.M()>500.)
            {
                m_histos.FillHisto(108+hoffset, jet2->calibratedEt(), weight, sysNum);
                foundOppositeEtaJetsM500 = true;
                break;
            }
        }
        if(foundOppositeEtaJetsM500) break;
    }
    if(!foundOppositeEtaJetsM500) m_histos.FillHisto(108+hoffset, 0., weight, sysNum);


    for(const auto& jet : m_sortedJets)
    {
        if(jet->calibratedEt()<20.) continue;
        m_histos.FillHisto(200+hoffset, jet->calibratedEt() , weight, sysNum);
        m_histos.FillHisto(201+hoffset, fabs(jet->eta()) , weight, sysNum);
    }
}


