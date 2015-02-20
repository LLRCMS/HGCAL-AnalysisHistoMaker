/**
 *  @file  AnalysisJets.cpp
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

#include "AnHiMaHGCAL/StudyJets/interface/AnalysisJets.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisJets::AnalysisJets():IAnalysis(),
    m_maxMassJets(0),
    m_vbfJets(0),
    m_random(561034)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisJets::~AnalysisJets()
/*****************************************************************/
{
    if(m_maxMassJets) delete m_maxMassJets;
}



/*****************************************************************/
bool AnalysisJets::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_egammaAlgo.initialize(event(), m_reader.params());
    m_jetAlgo.initialize(event(), m_reader.params());

    /// get BDT cuts
    TFile* fileBDTcuts = TFile::Open("/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/bdtCutsVsEta.root");
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
void AnalysisJets::execute()
/*****************************************************************/
{
    cerr<<"INFO: Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;

    m_emSeeds.clear();
    m_emClusters.clear();
    m_emClustersHard.clear();
    m_hadClusters.clear();
    m_hadClustersHard.clear();
    m_clusters.clear();
    m_clustersHard.clear();
    m_superClusters.clear();
    m_superClustersHard.clear();
    m_jets.clear();
    m_jetsHard.clear();

    m_maxMassJets = 0;
    m_vbfJets = 0;

    //-------------------------------------------
    // e/g algo
    // seeding
    m_egammaAlgo.seeding(event(), m_emSeeds);
    // clustering
    m_egammaAlgo.clustering(event(), m_emSeeds, m_emClusters);
    m_egammaAlgo.clustering(event(), m_emSeeds, m_emClustersHard, true, false);
    // cluster energy corrections
    m_egammaAlgo.clusterPileupSubtraction(m_emClusters);
    m_egammaAlgo.clusterThresholdCorrection(m_emClusters);
    //-------------------------------------------
    // Jet algo
    // clustering
    m_jetAlgo.hcalSeeding(event(), event().simhitsFH(), m_hadClusters);
    m_jetAlgo.hcalSeeding(event(), event().hardsimhitsFH(), m_hadClustersHard);
    //
    m_clusters = m_emClusters;
    for(const auto cl : m_hadClusters)
    {
        if(cl.calibratedEt()>1. && cl.nHits()>5) m_clusters.push_back(cl);
    }
    //m_clusters.insert(m_clusters.end(), m_hadClusters().begin(), m_hadClusters.end());
    m_clustersHard = m_emClustersHard;
    for(const auto cl : m_hadClustersHard)
    {
        if(cl.calibratedEt()>1. && cl.nHits()>5) m_clustersHard.push_back(cl);
    }
    //m_clustersHard.insert(m_clustersHard.end(), m_hadClustersHard().begin(), m_hadClustersHard.end());
    m_jetAlgo.superClustering(event(), m_clusters, m_superClusters, 0.4);
    m_jetAlgo.superClustering(event(), m_clustersHard, m_superClustersHard, 0.4);
    //
    vector<SimHit> simhits = event().simhits();
    simhits.insert(simhits.end(), event().simhitsFH().begin(), event().simhitsFH().end());
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets, 0.4);

    selectMaxMassJetPair();
    selectVBFJetPair();

    fillHistos();
}



/*****************************************************************/
void AnalysisJets::fillHistos()
/*****************************************************************/
{
    fillHistosAllJets();
    fillHistosVBFJets();
    fillHistosRandom();
    fillHistosTaus();
}

/*****************************************************************/
void AnalysisJets::fillHistosAllJets()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    // dijet
    if(m_maxMassJets)
    {
        double deta = fabs(m_maxMassJets->constituent(1).Eta() - m_maxMassJets->constituent(0).Eta());
        m_histos.FillHisto(100+hoffset, m_maxMassJets->M(), weight, sysNum);
        m_histos.FillHisto(101+hoffset, deta, weight, sysNum);
    }
    else
    {
        m_histos.FillHisto(100+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(101+hoffset, 0., weight, sysNum);
    }

    // hadronic clusters in event
    m_histos.FillHisto(21+hoffset, m_hadClusters.size(), weight, sysNum);
    m_histos.FillHisto(22+hoffset, m_hadClustersHard.size(), weight, sysNum);
    for(const auto& cluster : m_hadClusters)
    {
        m_histos.FillHisto(23+hoffset, cluster.calibratedEt(), weight, sysNum);
        m_histos.FillHisto(25+hoffset, cluster.calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(27+hoffset, cluster.nHits(), weight, sysNum);
        m_histos.FillHisto(29+hoffset, cluster.nLayers(4), weight, sysNum);
    }
    for(const auto& cluster : m_hadClustersHard)
    {
        m_histos.FillHisto(24+hoffset, cluster.calibratedEt(), weight, sysNum);
        m_histos.FillHisto(26+hoffset, cluster.calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(28+hoffset, cluster.nHits(), weight, sysNum);
        m_histos.FillHisto(30+hoffset, cluster.nLayers(4), weight, sysNum);
    }

    int nHADClustersAfterCuts = 0;
    for(const auto& cluster : m_hadClusters)
    {
        if(cluster.calibratedEt()>1. && cluster.nHits()>5) nHADClustersAfterCuts++;
    }
    m_histos.FillHisto(31+hoffset, nHADClustersAfterCuts, weight, sysNum);

    // Regions of interests
    m_histos.FillHisto(110+hoffset, m_superClusters.size(), weight, sysNum);

    // all jets
    m_histos.FillHisto(10+hoffset, event().genjets().size(), weight, sysNum);
    for(const auto& jet : event().genjets())
    {
        fillHistosJet(jet, sysNum, weight, hoffset);
    }


}


/*****************************************************************/
void AnalysisJets::fillHistosVBFJets()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 1000;


    if(!m_vbfJets) return;


    // vbf jet pair
    double deta = fabs(m_maxMassJets->constituent(1).Eta() - m_maxMassJets->constituent(0).Eta());
    m_histos.FillHisto(100+hoffset, m_maxMassJets->M(), weight, sysNum);
    m_histos.FillHisto(101+hoffset, deta, weight, sysNum);

    // vbf jets
    m_histos.FillHisto(10+hoffset, 2, weight, sysNum);
    for(unsigned i=0;i<2;i++)
    {
        const GenJet& jet = m_vbfJets->constituent(i);
        fillHistosJet(jet, sysNum, weight, hoffset);
    }
}

/*****************************************************************/
void AnalysisJets::fillHistosRandom()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 2000;

    int nDummyJets = 0;

    for(unsigned i=0;i<10;i++)
    {
        // choose random direction in HGCal
        double eta = m_random.Uniform(1.7,2.8);
        int side = m_random.Integer(2)*2-1;
        eta *= side;
        double phi = m_random.Uniform(-TMath::Pi(),TMath::Pi());

        double overlap = false;
        for(const auto& jet : event().genjets())
        {
            //if(jet.Pt()<15.) continue;
            double deta = eta - jet.Eta();
            double dphi =  TVector2::Phi_mpi_pi(phi - jet.Phi());
            double dr = sqrt(deta*deta+dphi*dphi);
            if(dr<0.4)
            {
                overlap = true;
                break;
            }
        }
        if(!overlap)
        {
            GenJet dummyJet;
            dummyJet.SetPtEtaPhiM(40., eta, phi, 0.);

            nDummyJets++;
            fillHistosJet(dummyJet, sysNum, weight, hoffset);
        }
    }
    m_histos.FillHisto(10+hoffset, nDummyJets, weight, sysNum);
}

/*****************************************************************/
void AnalysisJets::fillHistosTaus()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 3000;

    // dijet
    m_histos.FillHisto(100+hoffset, 0., weight, sysNum);
    m_histos.FillHisto(101+hoffset, 0., weight, sysNum);

    // all taus
    m_histos.FillHisto(10+hoffset, event().gentaus().size(), weight, sysNum);
    for(const auto& tau : event().gentaus())
    {
        fillHistosTau(tau, sysNum, weight, hoffset);
    }


}



/*****************************************************************/
void AnalysisJets::fillHistosJet(const GenJet& jet, short sysNum, float weight, int hoffset)
/*****************************************************************/
{
    // jet requirements
    if(fabs(jet.Eta())<1.5) return;
    if(fabs(jet.Eta())>3.) return;
    if(jet.Pt()<30.) return;


    m_histos.FillHisto(11+hoffset, jet.Pt()                     , weight, sysNum);
    m_histos.FillHisto(12+hoffset, fabs(jet.Eta())              , weight, sysNum);
    m_histos.FillHisto(13+hoffset, jet.E()                      , weight, sysNum);
    m_histos.FillHisto(14+hoffset, jet.emEnergy()               , weight, sysNum);
    m_histos.FillHisto(15+hoffset, jet.emEnergy()/jet.E()       , weight, sysNum);
    m_histos.FillHisto(16+hoffset, jet.hadEnergy()/jet.E()      , weight, sysNum);
    m_histos.FillHisto(17+hoffset, jet.invisibleEnergy()/jet.E(), weight, sysNum);
    m_histos.FillHisto(19+hoffset, 0                            , weight, sysNum);

    vector<const Tower*> clusters = matchClusters(jet); // match clusters to jet (including PU)
    // FIXME: decide if better without or with PU?
    vector<const Tower*> hadClusters = matchHadClusters(jet); // match had clusters to jet 
    vector<const SuperCluster*> rois = matchROI(jet); // match ROIs to jet 
    vector<const Tower*> jets = matchJet(jet); // match L1 jets to jet
    Tower coneEH = coneClustering(jet.Eta(), jet.Phi(), 0.4, true, true); // used to get reco jet eta,phi
    double detacone = jet.Eta()-coneEH.eta();
    double dphicone = TVector2::Phi_mpi_pi(jet.Phi()-coneEH.phi());
    double drcone = sqrt(detacone*detacone+dphicone*dphicone);

    if(clusters.size()>0)
    {
        m_histos.FillHisto(18+hoffset, jet.emEnergy()/jet.E()   , weight, sysNum);
        m_histos.FillHisto(20+hoffset, 0                        , weight, sysNum);
    }

    // electromagnetic cluster
    // sum of matched cluster energies
    double l1EmEnergy = 0.;
    for(const auto& cl : clusters) l1EmEnergy += cl->calibratedEnergy();

    m_histos.FillHisto(50+hoffset, clusters.size(), weight, sysNum);
    for(const auto& cl : clusters) 
    {
        m_histos.FillHisto(51+hoffset, cl->calibratedEnergy(), weight, sysNum);
    }

    m_histos.FillHisto(52+hoffset, (clusters.size()>0 ? clusters[0]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(53+hoffset, (clusters.size()>1 ? clusters[1]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(54+hoffset, (clusters.size()>2 ? clusters[2]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(55+hoffset, l1EmEnergy, weight, sysNum);


    // hadronic clusters
    double maxLayerEnergy = 0.;
    int maxLayer = 0;
    int firstLayer = 0;
    if(hadClusters.size()>0)
    {
        for(unsigned l=1;l<=12;l++)
        {
            double energy = hadClusters[0]->layerCalibratedEnergy(l, 4);
            if(energy>maxLayerEnergy)
            {
                maxLayer = l;
                maxLayerEnergy = energy;
            }
            if(firstLayer==0 && energy>0.) firstLayer = l;
        }
    }

    int nClusterAfterCuts = 0;
    int nClusterAfterEtCut = 0;
    for(const auto& cl : hadClusters)
    {
        if(cl->calibratedEt()>1.) nClusterAfterEtCut++;
        if(cl->calibratedEt()>1. && cl->nHits()>5) nClusterAfterCuts++;
    }

    //const Tower* leadingCluster = (hadClusters.size()>0 ? hadClusters[0] : 0);
    for(const auto& cl : hadClusters)
    {
        double dr = 999.;
        double drgen = 999.;
        double deta = cl->eta() - coneEH.eta();
        double dphi = TVector2::Phi_mpi_pi(cl->phi() - coneEH.phi());
        double detagen = cl->eta() - jet.Eta();
        double dphigen = TVector2::Phi_mpi_pi(cl->phi() - jet.Phi());
        dr = sqrt(deta*deta+dphi*dphi);
        drgen = sqrt(detagen*detagen+dphigen*dphigen);
        m_histos.FillHisto(64+hoffset, dr, weight*cl->calibratedEt(), sysNum);
        m_histos.FillHisto(65+hoffset, drgen, weight*cl->calibratedEt(), sysNum);
    }


    m_histos.FillHisto(56+hoffset, hadClusters.size(), weight, sysNum);
    m_histos.FillHisto(57+hoffset, (hadClusters.size()>0 ? hadClusters[0]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(58+hoffset, (hadClusters.size()>0 ? hadClusters[0]->calibratedEt() : 0.), weight, sysNum);
    m_histos.FillHisto(59+hoffset, (hadClusters.size()>0 ? hadClusters[0]->nHits() : 0), weight, sysNum);
    m_histos.FillHisto(60+hoffset, (hadClusters.size()>0 ? hadClusters[0]->nLayers(4) : 0), weight, sysNum);
    m_histos.FillHisto(61+hoffset, firstLayer, weight, sysNum);
    m_histos.FillHisto(62+hoffset, maxLayer, weight, sysNum);
    m_histos.FillHisto(63+hoffset, (hadClusters.size()>0 ? hadClusters[0]->longitudinalBarycenter(4) : 0.), weight, sysNum);
    m_histos.FillHisto(66+hoffset, drcone, weight, sysNum);
    m_histos.FillHisto(67+hoffset, nClusterAfterEtCut, weight, sysNum);
    m_histos.FillHisto(68+hoffset, nClusterAfterCuts, weight, sysNum);

    // compute cones in ECAL part. Around 0.4 E+H cone barycenter
    Tower cone005 = coneClustering(coneEH.eta(), coneEH.phi(), 0.05);
    Tower cone01  = coneClustering(coneEH.eta(), coneEH.phi(), 0.1);
    Tower cone02  = coneClustering(coneEH.eta(), coneEH.phi(), 0.2);
    Tower cone04  = coneClustering(coneEH.eta(), coneEH.phi(), 0.4);
    double l1EmEnergyCone005 = cone005.calibratedEnergy();
    double l1EmEnergyCone01  = cone01.calibratedEnergy();
    double l1EmEnergyCone02  = cone02.calibratedEnergy();
    double l1EmEnergyCone04  = cone04.calibratedEnergy();
    //
    Tower coneHard005 = coneClustering(coneEH.eta(), coneEH.phi(), 0.05, true);
    Tower coneHard01  = coneClustering(coneEH.eta(), coneEH.phi(), 0.1, true);
    Tower coneHard02  = coneClustering(coneEH.eta(), coneEH.phi(), 0.2, true);
    Tower coneHard04  = coneClustering(coneEH.eta(), coneEH.phi(), 0.4, true);
    double l1EmEnergyConeHard005 = coneHard005.calibratedEnergy();
    double l1EmEnergyConeHard01  = coneHard01.calibratedEnergy();
    double l1EmEnergyConeHard02  = coneHard02.calibratedEnergy();
    double l1EmEnergyConeHard04  = coneHard04.calibratedEnergy();

    m_histos.FillHisto(70+hoffset, l1EmEnergyCone005, weight, sysNum);
    m_histos.FillHisto(71+hoffset, l1EmEnergyCone01 , weight, sysNum);
    m_histos.FillHisto(72+hoffset, l1EmEnergyCone02 , weight, sysNum);
    m_histos.FillHisto(73+hoffset, l1EmEnergyCone04 , weight, sysNum);
    //
    m_histos.FillHisto(80+hoffset, l1EmEnergyConeHard005, weight, sysNum);
    m_histos.FillHisto(81+hoffset, l1EmEnergyConeHard01 , weight, sysNum);
    m_histos.FillHisto(82+hoffset, l1EmEnergyConeHard02 , weight, sysNum);
    m_histos.FillHisto(83+hoffset, l1EmEnergyConeHard04 , weight, sysNum);

    // Regions of interests
    m_histos.FillHisto(120+hoffset, rois.size(), weight, sysNum);
    for(const auto& roi : rois)
    {
        double dr = 999.;
        double drgen = 999.;
        double deta = roi->weightedEta() - coneEH.eta();
        double dphi = TVector2::Phi_mpi_pi(roi->weightedPhi() - coneEH.phi());
        double detagen = roi->weightedEta() - jet.Eta();
        double dphigen = TVector2::Phi_mpi_pi(roi->weightedPhi() - jet.Phi());
        dr = sqrt(deta*deta+dphi*dphi);
        drgen = sqrt(detagen*detagen+dphigen*dphigen);
        m_histos.FillHisto(121+hoffset, dr, weight*roi->et(), sysNum);
        m_histos.FillHisto(122+hoffset, drgen, weight*roi->et(), sysNum);
    }

    // L1 jets
    m_histos.FillHisto(150+hoffset, jets.size(), weight, sysNum);
    for(const auto& l1jet : jets)
    {
        double dr = 999.;
        double drgen = 999.;
        double deta = l1jet->eta() - coneEH.eta();
        double dphi = TVector2::Phi_mpi_pi(l1jet->phi() - coneEH.phi());
        double detagen = l1jet->eta() - jet.Eta();
        double dphigen = TVector2::Phi_mpi_pi(l1jet->phi() - jet.Phi());
        dr = sqrt(deta*deta+dphi*dphi);
        drgen = sqrt(detagen*detagen+dphigen*dphigen);
        m_histos.FillHisto(151+hoffset, dr, weight*l1jet->et(), sysNum);
        m_histos.FillHisto(152+hoffset, drgen, weight*l1jet->et(), sysNum);
    }
    if(rois.size()>0)
    {
        const SuperCluster* l1jet = rois[0];
        double minDR = 999.;
        for(unsigned i=0;i<rois.size();i++)
        {
            const SuperCluster* l1j = rois[i];
            double dr = 999.;
            double deta = l1j->weightedEta() - jet.Eta();
            double dphi = TVector2::Phi_mpi_pi(l1j->weightedPhi() - jet.Phi());
            dr = sqrt(deta*deta+dphi*dphi);
            if(dr<minDR)
            {
                minDR = dr;
                l1jet = l1j;
            }
        }
        //cout<<"minDR = "<<minDR<<"\n";
        vector<double> coneSizes = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4};
        double eta0 = l1jet->weightedEta();
        double phi0 = l1jet->weightedPhi();
        vector<Tower> l1jets = coneClustering(eta0, phi0, coneSizes, false);
        vector<Tower> l1jetsHard = coneClustering(eta0, phi0, coneSizes, true);
        m_histos.FillHisto(160+hoffset, l1jets[0].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(161+hoffset, l1jets[1].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(162+hoffset, l1jets[2].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(163+hoffset, l1jets[3].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(164+hoffset, l1jets[4].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(165+hoffset, l1jets[5].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(166+hoffset, l1jets[6].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(170+hoffset, l1jetsHard[0].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(171+hoffset, l1jetsHard[1].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(172+hoffset, l1jetsHard[2].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(173+hoffset, l1jetsHard[3].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(174+hoffset, l1jetsHard[4].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(175+hoffset, l1jetsHard[5].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(176+hoffset, l1jetsHard[6].calibratedEt(), weight, sysNum);
    }
    else
    {
        m_histos.FillHisto(160+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(161+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(162+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(163+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(164+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(165+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(166+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(170+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(171+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(172+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(173+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(174+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(175+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(176+hoffset, 0., weight, sysNum);
    }


    //
    // fill ntuple
    m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);
}

/*****************************************************************/
void AnalysisJets::fillHistosTau(const GenTau& tau, short sysNum, float weight, int hoffset)
/*****************************************************************/
{
    // tau requirements
    if(fabs(tau.Eta())<1.5) return;
    if(fabs(tau.Eta())>3.) return;
    if(tau.Pt()<20.) return;


    m_histos.FillHisto(11+hoffset, tau.Pt()                     , weight, sysNum);
    m_histos.FillHisto(12+hoffset, fabs(tau.Eta())              , weight, sysNum);
    m_histos.FillHisto(13+hoffset, tau.E()                      , weight, sysNum);
    m_histos.FillHisto(14+hoffset, 0.                           , weight, sysNum);
    m_histos.FillHisto(15+hoffset, 0.                           , weight, sysNum);
    m_histos.FillHisto(16+hoffset, 0.                           , weight, sysNum);
    m_histos.FillHisto(17+hoffset, 0.                           , weight, sysNum);
    m_histos.FillHisto(19+hoffset, tau.decay()                  , weight, sysNum);

    vector<const Tower*> clusters = matchClusters(tau); // match clusters to tau (including PU)
    vector<const Tower*> hadClusters = matchHadClusters(tau); // match had clusters to tau 
    vector<const SuperCluster*> rois = matchROI(tau); // match ROIs to tau
    vector<const Tower*> jets = matchJet(tau); // match L1 jets to jet
    Tower coneEH = coneClustering(tau.Eta(), tau.Phi(), 0.4, true, true); // used to get reco tau eta,phi
    double detacone = tau.Eta()-coneEH.eta();
    double dphicone = TVector2::Phi_mpi_pi(tau.Phi()-coneEH.phi());
    double drcone = sqrt(detacone*detacone+dphicone*dphicone);

    if(clusters.size()>0)
    {
        m_histos.FillHisto(18+hoffset, 0.         , weight, sysNum);
        m_histos.FillHisto(20+hoffset, tau.decay(), weight, sysNum);
    }

    // electromagnetic clusters
    // sum of matched cluster energies
    double l1EmEnergy = 0.;
    for(const auto& cl : clusters) l1EmEnergy += cl->calibratedEnergy();

    m_histos.FillHisto(50+hoffset, clusters.size(), weight, sysNum);
    for(const auto& cl : clusters) 
    {
        m_histos.FillHisto(51+hoffset, cl->calibratedEnergy(), weight, sysNum);
    }

    m_histos.FillHisto(52+hoffset, (clusters.size()>0 ? clusters[0]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(53+hoffset, (clusters.size()>1 ? clusters[1]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(54+hoffset, (clusters.size()>2 ? clusters[2]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(55+hoffset, l1EmEnergy, weight, sysNum);

    // hadronic clusters
    double maxLayerEnergy = 0.;
    int maxLayer = 0;
    int firstLayer = 0;
    if(hadClusters.size()>0)
    {
        for(unsigned l=1;l<=12;l++)
        {
            double energy = hadClusters[0]->layerCalibratedEnergy(l, 4);
            if(energy>maxLayerEnergy)
            {
                maxLayer = l;
                maxLayerEnergy = energy;
            }
            if(firstLayer==0 && energy>0.) firstLayer = l;
        }
    }

    int nClusterAfterCuts = 0;
    int nClusterAfterEtCut = 0;
    for(const auto& cl : hadClusters)
    {
        if(cl->calibratedEt()>1.) nClusterAfterEtCut++;
        if(cl->calibratedEt()>1. && cl->nHits()>5) nClusterAfterCuts++;
    }

    //const Tower* leadingCluster = (hadClusters.size()>0 ? hadClusters[0] : 0);
    for(const auto& cl : hadClusters)
    {
        double dr = 999.;
        double drgen = 999.;
        double deta = cl->eta() - coneEH.eta();
        double dphi = TVector2::Phi_mpi_pi(cl->phi() - coneEH.phi());
        double detagen = cl->eta() - tau.Eta();
        double dphigen = TVector2::Phi_mpi_pi(cl->phi() - tau.Phi());
        dr = sqrt(deta*deta+dphi*dphi);
        drgen = sqrt(detagen*detagen+dphigen*dphigen);
        m_histos.FillHisto(64+hoffset, dr, weight*cl->calibratedEt(), sysNum);
        m_histos.FillHisto(65+hoffset, drgen, weight*cl->calibratedEt(), sysNum);
    }


    m_histos.FillHisto(56+hoffset, hadClusters.size(), weight, sysNum);
    m_histos.FillHisto(57+hoffset, (hadClusters.size()>0 ? hadClusters[0]->calibratedEnergy() : 0.), weight, sysNum);
    m_histos.FillHisto(58+hoffset, (hadClusters.size()>0 ? hadClusters[0]->calibratedEt() : 0.), weight, sysNum);
    m_histos.FillHisto(59+hoffset, (hadClusters.size()>0 ? hadClusters[0]->nHits() : 0), weight, sysNum);
    m_histos.FillHisto(60+hoffset, (hadClusters.size()>0 ? hadClusters[0]->nLayers(4) : 0), weight, sysNum);
    m_histos.FillHisto(61+hoffset, firstLayer, weight, sysNum);
    m_histos.FillHisto(62+hoffset, maxLayer, weight, sysNum);
    m_histos.FillHisto(63+hoffset, (hadClusters.size()>0 ? hadClusters[0]->longitudinalBarycenter(4) : 0.), weight, sysNum);

    m_histos.FillHisto(66+hoffset, drcone, weight, sysNum);
    m_histos.FillHisto(67+hoffset, nClusterAfterEtCut, weight, sysNum);
    m_histos.FillHisto(68+hoffset, nClusterAfterCuts, weight, sysNum);


    // compute cones in ECAL part. Around 0.4 E+H cone barycenter
    Tower cone005 = coneClustering(coneEH.eta(), coneEH.phi(), 0.05);
    Tower cone01  = coneClustering(coneEH.eta(), coneEH.phi(), 0.1);
    Tower cone02  = coneClustering(coneEH.eta(), coneEH.phi(), 0.2);
    Tower cone04  = coneClustering(coneEH.eta(), coneEH.phi(), 0.4);
    double l1EmEnergyCone005 = cone005.calibratedEnergy();
    double l1EmEnergyCone01  = cone01.calibratedEnergy();
    double l1EmEnergyCone02  = cone02.calibratedEnergy();
    double l1EmEnergyCone04  = cone04.calibratedEnergy();
    //
    Tower coneHard005 = coneClustering(coneEH.eta(), coneEH.phi(), 0.05, true);
    Tower coneHard01  = coneClustering(coneEH.eta(), coneEH.phi(), 0.1, true);
    Tower coneHard02  = coneClustering(coneEH.eta(), coneEH.phi(), 0.2, true);
    Tower coneHard04  = coneClustering(coneEH.eta(), coneEH.phi(), 0.4, true);
    double l1EmEnergyConeHard005 = coneHard005.calibratedEnergy();
    double l1EmEnergyConeHard01  = coneHard01.calibratedEnergy();
    double l1EmEnergyConeHard02  = coneHard02.calibratedEnergy();
    double l1EmEnergyConeHard04  = coneHard04.calibratedEnergy();

    m_histos.FillHisto(70+hoffset, l1EmEnergyCone005, weight, sysNum);
    m_histos.FillHisto(71+hoffset, l1EmEnergyCone01 , weight, sysNum);
    m_histos.FillHisto(72+hoffset, l1EmEnergyCone02 , weight, sysNum);
    m_histos.FillHisto(73+hoffset, l1EmEnergyCone04 , weight, sysNum);
    //
    m_histos.FillHisto(80+hoffset, l1EmEnergyConeHard005, weight, sysNum);
    m_histos.FillHisto(81+hoffset, l1EmEnergyConeHard01 , weight, sysNum);
    m_histos.FillHisto(82+hoffset, l1EmEnergyConeHard02 , weight, sysNum);
    m_histos.FillHisto(83+hoffset, l1EmEnergyConeHard04 , weight, sysNum);

    // Regions of interests
    m_histos.FillHisto(120+hoffset, rois.size(), weight, sysNum);
    for(const auto& roi : rois)
    {
        double dr = 999.;
        double drgen = 999.;
        double deta = roi->weightedEta() - coneEH.eta();
        double dphi = TVector2::Phi_mpi_pi(roi->weightedPhi() - coneEH.phi());
        double detagen = roi->weightedEta() - tau.Eta();
        double dphigen = TVector2::Phi_mpi_pi(roi->weightedPhi() - tau.Phi());
        dr = sqrt(deta*deta+dphi*dphi);
        drgen = sqrt(detagen*detagen+dphigen*dphigen);
        m_histos.FillHisto(121+hoffset, dr, weight*roi->et(), sysNum);
        m_histos.FillHisto(122+hoffset, drgen, weight*roi->et(), sysNum);
    }

    // L1 jets
    m_histos.FillHisto(150+hoffset, jets.size(), weight, sysNum);
    for(const auto& l1jet : jets)
    {
        double dr = 999.;
        double drgen = 999.;
        double deta = l1jet->eta() - coneEH.eta();
        double dphi = TVector2::Phi_mpi_pi(l1jet->phi() - coneEH.phi());
        double detagen = l1jet->eta() - tau.Eta();
        double dphigen = TVector2::Phi_mpi_pi(l1jet->phi() - tau.Phi());
        dr = sqrt(deta*deta+dphi*dphi);
        drgen = sqrt(detagen*detagen+dphigen*dphigen);
        m_histos.FillHisto(151+hoffset, dr, weight*l1jet->et(), sysNum);
        m_histos.FillHisto(152+hoffset, drgen, weight*l1jet->et(), sysNum);
    }
    if(rois.size()>0)
    {
        const SuperCluster* l1jet = rois[0];
        double minDR = 999.;
        for(unsigned i=0;i<rois.size();i++)
        {
            const SuperCluster* l1j = rois[i];
            double dr = 999.;
            double deta = l1j->weightedEta() - tau.Eta();
            double dphi = TVector2::Phi_mpi_pi(l1j->weightedPhi() - tau.Phi());
            dr = sqrt(deta*deta+dphi*dphi);
            if(dr<minDR)
            {
                minDR = dr;
                l1jet = l1j;
            }
        }
        vector<double> coneSizes = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4};
        double eta0 = l1jet->weightedEta();
        double phi0 = l1jet->weightedPhi();
        vector<Tower> l1jets = coneClustering(eta0, phi0, coneSizes, false);
        vector<Tower> l1jetsHard = coneClustering(eta0, phi0, coneSizes, true);
        m_histos.FillHisto(160+hoffset, l1jets[0].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(161+hoffset, l1jets[1].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(162+hoffset, l1jets[2].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(163+hoffset, l1jets[3].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(164+hoffset, l1jets[4].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(165+hoffset, l1jets[5].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(166+hoffset, l1jets[6].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(170+hoffset, l1jetsHard[0].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(171+hoffset, l1jetsHard[1].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(172+hoffset, l1jetsHard[2].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(173+hoffset, l1jetsHard[3].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(174+hoffset, l1jetsHard[4].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(175+hoffset, l1jetsHard[5].calibratedEt(), weight, sysNum);
        m_histos.FillHisto(176+hoffset, l1jetsHard[6].calibratedEt(), weight, sysNum);
    }
    else
    {
        m_histos.FillHisto(160+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(161+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(162+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(163+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(164+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(165+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(166+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(170+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(171+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(172+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(173+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(174+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(175+hoffset, 0., weight, sysNum);
        m_histos.FillHisto(176+hoffset, 0., weight, sysNum);
    }

    //
    // fill ntuple
    m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);
}


/*****************************************************************/
void AnalysisJets::selectMaxMassJetPair()
/*****************************************************************/
{
    if(m_maxMassJets) delete m_maxMassJets;
    m_maxMassJets = 0;

    vector< TParticlePair<GenJet> > jetPairs;
    unsigned nJets = event().genjets().size();
    for(int i=0;i<(int)nJets-1;i++)
    {
        const GenJet& jet1 = event().genjets()[i];
        if(jet1.Pt()<30.) continue;
        if(overlapWithLeptonOrNeutrino(jet1)) continue;
        for(int j=i+1;j<(int)nJets;j++)
        {
            const GenJet& jet2 = event().genjets()[j];
            if(jet2.Pt()<30.) continue;
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
void AnalysisJets::selectVBFJetPair()
/*****************************************************************/
{
    m_vbfJets = 0;
    if(!m_maxMassJets) return;
    if(m_maxMassJets->M()<250.) return;
    if(fabs(m_maxMassJets->constituent(1).Eta() - m_maxMassJets->constituent(0).Eta())<2.5) return; 

    m_vbfJets = m_maxMassJets;
}


/*****************************************************************/
bool AnalysisJets::overlapWithLeptonOrNeutrino(const GenJet& jet)
/*****************************************************************/
{
    bool overlap = false;
    //int pid = 0;
    //double pt = 0.;
    for(const auto& p : event().genparticles())
    {
        if(abs(p.id())!=15) continue; // select taus
        //if(p.status()!=1) continue;
        //if(abs(p.id())!=12 && 
        //   abs(p.id())!=13 &&
        //   abs(p.id())!=14 &&
        //   abs(p.id())!=15 &&
        //   abs(p.id())!=16
        //  ) continue;

        if(p.DeltaR(jet)<0.4)
        {
            overlap = true;
            //pid = p.id();
            //pt = p.Pt();
        }
    }

    return overlap;
}


/*****************************************************************/
vector<const Tower*> AnalysisJets::matchClusters(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> clusters;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // with PU
    for(const auto& cl : m_emClusters )
    {
        double deta = (cl.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) clusters.push_back(&cl);
    }
    sort(clusters.begin(),clusters.end(), AnHiMa::towerSortE);
    return clusters;
}

/*****************************************************************/
vector<const Tower*> AnalysisJets::matchClustersHard(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> clustersHard;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // w/o PU
    for(const auto& cl : m_emClustersHard )
    {
        double deta = (cl.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) clustersHard.push_back(&cl);
    }
    sort(clustersHard.begin(),clustersHard.end(), AnHiMa::towerSortE);
    return clustersHard;
}


/*****************************************************************/
vector<const Tower*> AnalysisJets::matchHadClusters(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> clusters;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // with PU
    for(const auto& cl : m_hadClusters )
    {
        double deta = (cl.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) clusters.push_back(&cl);
    }
    sort(clusters.begin(),clusters.end(), AnHiMa::towerSortE);
    return clusters;
}

/*****************************************************************/
vector<const Tower*> AnalysisJets::matchHadClustersHard(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> clustersHard;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // w/o PU
    for(const auto& cl : m_hadClustersHard )
    {
        double deta = (cl.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) clustersHard.push_back(&cl);
    }
    sort(clustersHard.begin(),clustersHard.end(), AnHiMa::towerSortE);
    return clustersHard;
}

/*****************************************************************/
vector<const SuperCluster*> AnalysisJets::matchROI(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const SuperCluster*> rois;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // with PU
    for(const auto& cl : m_superClusters)
    {
        double deta = (cl.weightedEta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.weightedPhi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) rois.push_back(&cl);
    }
    sort(rois.begin(),rois.end(), AnHiMa::superClusterSort);
    return rois;
}

/*****************************************************************/
vector<const SuperCluster*> AnalysisJets::matchROIHard(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const SuperCluster*> rois;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // without PU
    for(const auto& cl : m_superClustersHard)
    {
        double deta = (cl.weightedEta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.weightedPhi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) rois.push_back(&cl);
    }
    sort(rois.begin(),rois.end(), AnHiMa::superClusterSort);
    return rois;
}

/*****************************************************************/
vector<const Tower*> AnalysisJets::matchJet(const TLorentzVector& jet)
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
        if(dr<0.4) jets.push_back(&cl);
    }
    sort(jets.begin(),jets.end(), AnHiMa::towerSort);
    return jets;
}

/*****************************************************************/
vector<const Tower*> AnalysisJets::matchJetHard(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> jets;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // with PU
    for(const auto& cl : m_jetsHard)
    {
        double deta = (cl.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(cl.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.4) jets.push_back(&cl);
    }
    sort(jets.begin(),jets.end(), AnHiMa::towerSort);
    return jets;
}


/*****************************************************************/
Tower AnalysisJets::coneClustering(double eta0, double phi0, double coneSize, bool EH, bool hard)
/*****************************************************************/
{
    vector<SimHit> simhits = (hard ? event().hardsimhits() : event().simhits());
    if(EH)
    {
        const vector<SimHit>& simhitsFH = (hard ? event().hardsimhitsFH() : event().simhitsFH());
        simhits.insert(simhits.end(), simhitsFH.begin(), simhitsFH.end());
    }

    // to calibrate clusters
    TowerCalibrator towerCalibrator;

    // Build clusters around specified direction
    Tower cluster;
    for(const auto& hit: simhits)
    {
        float eta = hit.eta();
        float phi = hit.phi();
        double deta = (eta-eta0);
        double dphi =  TVector2::Phi_mpi_pi(phi-phi0);
        float dr = sqrt(deta*deta+dphi*dphi);
        if(dr>coneSize) continue;
        cluster.addHit(hit);
    }
    towerCalibrator.calibrate(cluster);
    return cluster;
}


/*****************************************************************/
vector<Tower> AnalysisJets::coneClustering(double eta0, double phi0, const vector<double>& coneSizes, bool hard)
/*****************************************************************/
{
    vector<SimHit> simhits = (hard ? event().hardsimhits() : event().simhits());
    const vector<SimHit>& simhitsFH = (hard ? event().hardsimhitsFH() : event().simhitsFH());
    simhits.insert(simhits.end(), simhitsFH.begin(), simhitsFH.end());

    // to calibrate clusters
    TowerCalibrator towerCalibrator;

    // Build clusters around specified direction
    vector<Tower> clusters(coneSizes.size());
    for(const auto& hit: simhits)
    {
        float eta = hit.eta();
        float phi = hit.phi();
        double deta = (eta-eta0);
        double dphi =  TVector2::Phi_mpi_pi(phi-phi0);
        float dr = sqrt(deta*deta+dphi*dphi);
        for(unsigned i=0;i<coneSizes.size();i++)
        {
            if(dr>coneSizes[i]) continue;
            clusters[i].addHit(hit);
        }
    }
    for(auto& cluster : clusters) towerCalibrator.calibrate(cluster);
    return clusters;
}

