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
    m_matchedSuperCluster(0),
    m_matchedIdealCluster(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisEfficiency::~AnalysisEfficiency()
/*****************************************************************/
{
    for(int i=0;i<1;++i)
    {
         m_failAndPassTrees[i]->Write();
    }

}



/*****************************************************************/
bool AnalysisEfficiency::initialize(const string& parameterFile)
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


    ///// turn-on trees
    m_thresholds = {25,27,30,35,40};
    m_failAndPassTrees = new TTree*[1];
    m_failAndPassTrees[0] = new TTree("TurnOn_L1SuperCluster_failAndPassTree","TurnOn_L1SuperCluster_failAndPassTree");
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
            stringstream passID, passIDtype;
            passID << "l1_passID_" << th;
            passIDtype << "l1_passID_" << th<<"/I";
            m_failAndPassIDBits[i][th] = 0;
            tree->Branch(passID.str().c_str(), &m_failAndPassIDBits[i][th], passIDtype.str().c_str());
        }
    }


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
            // ideal clustering
            m_idealClusters.clear();
            m_egammaAlgo.idealClustering(event(), m_idealClusters, m_genParticle->Eta(), m_genParticle->Phi());
            m_matchedIdealCluster = (m_idealClusters.size()>0 ? &m_idealClusters[0] : 0);
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
        m_histos.FillHisto(13+hoffset, m_genParticle->Pt(), weight, sysNum);
    }
    //
    if(m_matchedSeed)
    {
        m_histos.FillHisto(20+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(21+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(22+hoffset, event().npu(), weight, sysNum);
        m_histos.FillHisto(23+hoffset, m_genParticle->Pt(), weight, sysNum);
    }
    //
    if(m_matchedCluster && m_matchedCluster->calibratedEt()>5.)
    {
        m_histos.FillHisto(30+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(31+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(32+hoffset, event().npu(), weight, sysNum);
        m_histos.FillHisto(33+hoffset, m_genParticle->Pt(), weight, sysNum);
    }
    //
    if(m_matchedSuperCluster && m_matchedSuperCluster->et()>5.)
    {
        m_histos.FillHisto(40+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(41+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(42+hoffset, event().npu(), weight, sysNum);
        m_histos.FillHisto(43+hoffset, m_genParticle->Pt(), weight, sysNum);

        string bdtName = "BDTG_QCD";
        if(m_egammaAlgo.useHalfLayers()) bdtName = "BDTG_QCD_HalfLayers";
        double bdt = m_egammaAlgo.bdtOutput(*m_matchedSuperCluster->cluster(0), "BDTG_QCD");
        bool pass995 = (bdt>=m_bdtCuts[995]->Eval(fabs(m_matchedSuperCluster->eta())));
        bool pass99  = (bdt>=m_bdtCuts[99] ->Eval(fabs(m_matchedSuperCluster->eta())));
        bool pass985 = (bdt>=m_bdtCuts[985]->Eval(fabs(m_matchedSuperCluster->eta())));
        bool pass98  = (bdt>=m_bdtCuts[98] ->Eval(fabs(m_matchedSuperCluster->eta())));
        bool pass975 = (bdt>=m_bdtCuts[975]->Eval(fabs(m_matchedSuperCluster->eta())));
        bool pass97  = (bdt>=m_bdtCuts[97] ->Eval(fabs(m_matchedSuperCluster->eta())));
        //double reducedPhiSC = (m_matchedSuperCluster->phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
        //double localPhiSC = (reducedPhiSC - floor(reducedPhiSC))*(2*TMath::Pi()/18.);
        //bool regionLoose = (fabs(m_matchedSuperCluster->eta())<2.2 || m_matchedSuperCluster->et()<30. || localPhiSC<0.05 || localPhiSC>(2*TMath::Pi()/18.-0.05));
        if(pass995)
        {
            m_histos.FillHisto(60+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(61+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(62+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(63+hoffset, m_genParticle->Pt(), weight, sysNum);
        }
        if(pass99)
        {
            m_histos.FillHisto(70+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(71+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(72+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(73+hoffset, m_genParticle->Pt(), weight, sysNum);
        }
        if(pass985)
        {
            m_histos.FillHisto(80+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(81+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(82+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(83+hoffset, m_genParticle->Pt(), weight, sysNum);
        }
        if(pass98)
        {
            m_histos.FillHisto(90+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(91+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(92+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(93+hoffset, m_genParticle->Pt(), weight, sysNum);
        }
        if(pass975)
        {
            m_histos.FillHisto(100+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(101+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(102+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(103+hoffset, m_genParticle->Pt(), weight, sysNum);
        }
        if(pass97)
        {
            m_histos.FillHisto(110+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(111+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(112+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(113+hoffset, m_genParticle->Pt(), weight, sysNum);
        }
        if( (m_matchedSuperCluster->et()<30. && pass99) || (m_matchedSuperCluster->et()>30. && pass98))
        {
            m_histos.FillHisto(120+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
            m_histos.FillHisto(121+hoffset, localPhi, weight, sysNum);
            m_histos.FillHisto(122+hoffset, event().npu(), weight, sysNum);
            m_histos.FillHisto(123+hoffset, m_genParticle->Pt(), weight, sysNum);
        }


        // fill turn-on tree
        m_gen_pt[0] = m_genParticle->Pt();
        m_gen_eta[0] = m_genParticle->Eta();
        bool passID = ( (m_matchedSuperCluster->et()<30. && pass99) || (m_matchedSuperCluster->et()>30. && pass98) );
        for(unsigned i=0;i<m_thresholds.size();i++)
        {
            int th = m_thresholds[i];
            m_failAndPassBits[0][th] = (m_matchedSuperCluster->et()>=(double)th);
            m_failAndPassIDBits[0][th] = (m_matchedSuperCluster->et()>=(double)th && passID );
        }


        m_failAndPassTrees[0]->Fill();
    }
    //
    if(m_matchedIdealCluster && m_matchedIdealCluster->calibratedEt()>10.)
    {
        m_histos.FillHisto(50+hoffset, fabs(m_genParticle->Eta()), weight, sysNum);
        m_histos.FillHisto(51+hoffset, localPhi, weight, sysNum);
        m_histos.FillHisto(52+hoffset, event().npu(), weight, sysNum);
        m_histos.FillHisto(53+hoffset, m_genParticle->Pt(), weight, sysNum);
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







