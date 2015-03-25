/**
 *  @file  AnalysisJetEfficiency.cpp
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

#include "AnHiMaHGCAL/TriggerJetEfficiency/interface/AnalysisJetEfficiency.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisJetEfficiency::AnalysisJetEfficiency():IAnalysis(),
    m_maxMassJets(0),
    m_vbfJets(0),
    m_genJet(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisJetEfficiency::~AnalysisJetEfficiency()
/*****************************************************************/
{
    if(m_maxMassJets) delete m_maxMassJets;

    for(int i=0;i<1;++i)
    {
         m_failAndPassTrees[i]->Write();
    }

}



/*****************************************************************/
bool AnalysisJetEfficiency::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // will built projected hits onto the HGCEE layer 1
    string projectionMappingFile(m_reader.params().GetValue("ProjectionMapping", ""));
    event().initializeHitProjection(projectionMappingFile);


    // Read parameters
    //string pileupParamsFile = m_reader.params().GetValue("PileupParams", "/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC19/src/AnHiMaHGCAL/ElectronClusterThreshold/data/thresholdParameters.txt");

    m_egammaAlgo.initialize(event(), m_reader.params());
    m_jetAlgo.initialize(event(), m_reader.params());

    ///// turn-on trees
    m_thresholds = {30, 40, 60, 80, 100, 150};
    m_failAndPassTrees = new TTree*[1];
    m_failAndPassTrees[0] = new TTree("TurnOn_L1Jet_failAndPassTree","TurnOn_L1Jet_failAndPassTree");
    m_outputFile->cd();
    for(int i=0;i<1;++i)
    {
        //cerr<<i<<"\n";
        TTree* tree = m_failAndPassTrees[i];
        m_gen_pt[i] = 0.;
        m_gen_eta[i] = 0.;
        stringstream genptvar, genptvartype, genetavar, genetavartype;
        genptvar << "gen_pt";
        genptvartype << "gen_pt/D";
        genetavar << "gen_eta";
        genetavartype << "gen_eta/D";
        tree->Branch(genptvar.str().c_str(), &m_gen_pt[i], genptvartype.str().c_str());
        tree->Branch(genetavar.str().c_str(), &m_gen_eta[i], genetavartype.str().c_str());
        for(unsigned j=0;j<m_thresholds.size();j++)
        {
            int th = m_thresholds[j];
            stringstream pass, passtype;
            pass << "l1_pass_" << th;
            passtype << "l1_pass_" << th<<"/I";
            m_failAndPassBits[i][th] = 0;
            tree->Branch(pass.str().c_str(), &m_failAndPassBits[i][th], passtype.str().c_str());
            //
            stringstream passID, passIDtype;
            passID << "l1_passID_" << th;
            passIDtype << "l1_passID_" << th<<"/I";
            m_failAndPassIDBits[i][th] = 0;
        }
    }


    return true;
}


/*****************************************************************/
void AnalysisJetEfficiency::execute()
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
    m_jets.clear();
    m_jetCores.clear();

    m_maxMassJets = 0;
    m_vbfJets = 0;

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
    for(const auto cl : m_hadClusters)
    {
        if(cl.calibratedEt()>1. && cl.nHits()>5) m_clusters.push_back(cl);
    }
    m_jetAlgo.superClustering(event(), m_clusters, m_superClusters, 0.4);
    //
    vector<SimHit> simhits = event().simhits();
    simhits.insert(simhits.end(), event().simhitsFH().begin(), event().simhitsFH().end());
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets, 0.2); // standard cone
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jetCores, 0.1); // core

    m_jetAlgo.pileupSubtraction(m_jets, "0.2");
    m_jetAlgo.pileupSubtraction(m_jetCores, "0.1");
    m_jetAlgo.thresholdCorrection(m_jets, "0.2");
    m_jetAlgo.thresholdCorrection(m_jetCores, "0.1");
    m_jetAlgo.truthCorrection(m_jets, m_jetCores);

    selectMaxMassJetPair();
    selectVBFJetPair();

    if(m_vbfJets)
    {
        for(unsigned i=0;i<2;i++)
        {
            const GenJet& jet = (*m_vbfJets)[i];

            // jet requirements
            if(fabs(jet.Eta())<1.6) continue;
            if(fabs(jet.Eta())>2.9) continue;
            if(jet.Pt()<15.) continue;

            m_genJet = &(*m_vbfJets)[i];
            fillHistos();

        }
    }

}

/*****************************************************************/
void AnalysisJetEfficiency::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

    vector<const Tower*> l1jets = matchJet(*m_genJet);
    const Tower* l1jet          = (l1jets.size()>0 ? l1jets[0] : 0);

    double reducedPhi = (m_genJet->Phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
    double localPhi = (reducedPhi - floor(reducedPhi))*(2*TMath::Pi()/18.);

    m_histos.FillHisto(10+hoffset, fabs(m_genJet->Eta()), weight, sysNum);
    m_histos.FillHisto(11+hoffset, localPhi, weight, sysNum);
    m_histos.FillHisto(12+hoffset, event().npu(), weight, sysNum);
    m_histos.FillHisto(13+hoffset, m_genJet->Et(), weight, sysNum);
    if(l1jet && l1jet->calibratedEt()>5.)
    {
        m_histos.FillHisto(20+hoffset, fabs(m_genJet->Eta()), weight, sysNum);
        m_histos.FillHisto(21+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(22+hoffset, event().npu(), weight, sysNum);
        m_histos.FillHisto(23+hoffset, m_genJet->Et(), weight, sysNum);

        // fill turn-on tree
        m_gen_pt[0] = m_genJet->Et();
        m_gen_eta[0] = m_genJet->Eta();
        bool passID = true;
        for(unsigned i=0;i<m_thresholds.size();i++)
        {
            int th = m_thresholds[i];
            m_failAndPassBits[0][th]   = (l1jet->calibratedEt()>=(double)th);
            m_failAndPassIDBits[0][th] = (l1jet->calibratedEt()>=(double)th && passID);
        }

    }
    m_failAndPassTrees[0]->Fill();

}


/*****************************************************************/
void AnalysisJetEfficiency::selectMaxMassJetPair()
/*****************************************************************/
{
    if(m_maxMassJets) delete m_maxMassJets;
    m_maxMassJets = 0;

    vector< TParticlePair<GenJet> > jetPairs;
    unsigned nJets = event().genjets().size();
    for(int i=0;i<(int)nJets-1;i++)
    {
        const GenJet& jet1 = event().genjets()[i];
        if(jet1.Pt()<15.) continue;
        if(overlapWithLeptonOrNeutrino(jet1)) continue;
        for(int j=i+1;j<(int)nJets;j++)
        {
            const GenJet& jet2 = event().genjets()[j];
            if(jet2.Pt()<15.) continue;
            if(overlapWithLeptonOrNeutrino(jet2)) continue;
            jetPairs.push_back( TParticlePair<GenJet>(jet1,jet2) );
        }
    }

    const TParticlePair<GenJet>* maxMassPair = 0;
    for(const auto& jetPair : jetPairs)
    {
        if(!maxMassPair || jetPair.M()>maxMassPair->M())
        {
            maxMassPair = &jetPair;
        }
    }

    if(maxMassPair) m_maxMassJets = new TParticlePair<GenJet>(*maxMassPair);
}



/*****************************************************************/
void AnalysisJetEfficiency::selectVBFJetPair()
/*****************************************************************/
{
    m_vbfJets = 0;
    if(!m_maxMassJets) return;
    if(m_maxMassJets->M()<250.) return;
    if(m_maxMassJets->constituent(1).Eta()*m_maxMassJets->constituent(0).Eta()>0.) return; 
    if(fabs(m_maxMassJets->constituent(1).Eta() - m_maxMassJets->constituent(0).Eta())<2.5) return; 

    m_vbfJets = m_maxMassJets;
}


/*****************************************************************/
bool AnalysisJetEfficiency::overlapWithLeptonOrNeutrino(const GenJet& jet)
/*****************************************************************/
{
    bool overlap = false;
    for(const auto& p : event().genparticles())
    {
        if(abs(p.id())!=15) continue; // select taus

        if(p.DeltaR(jet)<0.4)
        {
            overlap = true;
        }
    }

    return overlap;
}



/*****************************************************************/
vector<const Tower*> AnalysisJetEfficiency::matchJet(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> jets;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // with PU
    for(const auto& cl : m_jets)
    {
        double deta = (cl.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2) jets.push_back(&cl);
    }
    sort(jets.begin(),jets.end(), AnHiMa::towerSort);
    return jets;
}




/*****************************************************************/
const Tower* AnalysisJetEfficiency::matchJetCore(const Tower& l1jet)
/*****************************************************************/
{
    const Tower* core = 0;
    double drMin = 99999.;
    double eta = l1jet.eta();
    double phi = l1jet.phi();
    for(const auto& c : m_jetCores)
    {
        double deta = (c.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(c.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.05 && dr<drMin)
        {
            if(core) cout<<"WARNING: L1 jet matched to several jet cores\n";
            core = &c;
            drMin = dr;
        }
    }
    if(!core)
    {
        cout<<"WARNING: cannot find core jet\n";
    }
    return core;
}










