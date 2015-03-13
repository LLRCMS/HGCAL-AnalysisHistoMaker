/**
 *  @file  AnalysisJetCalibration.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    02/02/2015
 *
 *  @internal
 *     Created :  02/02/2015
 * Last update :  02/02/2015 11:24:52
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */






#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/TriggerJetCalibration/interface/AnalysisJetCalibration.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisJetCalibration::AnalysisJetCalibration():IAnalysis(),
    m_maxMassJets(0),
    m_vbfJets(0),
    m_random(561034)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisJetCalibration::~AnalysisJetCalibration()
/*****************************************************************/
{
    if(m_maxMassJets) delete m_maxMassJets;
}



/*****************************************************************/
bool AnalysisJetCalibration::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_egammaAlgo.initialize(event(), m_reader.params());
    m_jetAlgo.initialize(event(), m_reader.params());


    return true;
}


/*****************************************************************/
void AnalysisJetCalibration::execute()
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
    m_jetCores.clear();
    m_jets04.clear();
    m_jetsHard.clear();
    m_jetCoresHard.clear();
    m_jets04Hard.clear();
    m_jetsHardNoThreshold.clear();
    m_jetCoresHardNoThreshold.clear();
    m_jets04HardNoThreshold.clear();

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
    vector<SimHit> simhitshard = event().hardsimhits();
    simhitshard.insert(simhitshard.end(), event().hardsimhitsFH().begin(), event().hardsimhitsFH().end());
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets, 0.2); // standard cone
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jetCores, 0.1); // core
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets04, 0.4); // 0.4
    m_jetAlgo.coneClustering(event(), simhitshard, m_superClusters, m_jetsHard, 0.2); // cone on hard hits
    m_jetAlgo.coneClustering(event(), simhitshard, m_superClusters, m_jetCoresHard, 0.1); // core on hard hits
    m_jetAlgo.coneClustering(event(), simhitshard, m_superClusters, m_jets04Hard, 0.4); // 0.4 on hard hits
    m_jetAlgo.coneClustering(event(), simhitshard, m_superClusters, m_jetsHardNoThreshold, 0.2, false); // cone on hard hits w/o threshold
    m_jetAlgo.coneClustering(event(), simhitshard, m_superClusters, m_jetCoresHardNoThreshold, 0.1, false); // core on hard hits w/o threshold
    m_jetAlgo.coneClustering(event(), simhitshard, m_superClusters, m_jets04HardNoThreshold, 0.4, false); // 0.4 on hard hits w/o threshold

    selectMaxMassJetPair();
    selectVBFJetPair();

    fillHistos();
}



/*****************************************************************/
void AnalysisJetCalibration::fillHistos()
/*****************************************************************/
{
    fillHistosVBFJets(0);
    m_jetAlgo.pileupSubtraction(m_jets, "0.2");
    m_jetAlgo.pileupSubtraction(m_jetCores, "0.1");
    m_jetAlgo.pileupSubtraction(m_jets04, "0.4");
    fillHistosVBFJets(1000);
    m_jetAlgo.thresholdCorrection(m_jets, "0.2");
    m_jetAlgo.thresholdCorrection(m_jetCores, "0.1");
    m_jetAlgo.thresholdCorrection(m_jets04, "0.4");
    fillHistosVBFJets(2000);
    m_jetAlgo.truthCorrection(m_jets, m_jetCores);
    fillHistosVBFJets(3000);
}


/*****************************************************************/
void AnalysisJetCalibration::fillHistosVBFJets(int hoffset)
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;



    if(!m_vbfJets) return;
    for(unsigned i=0;i<2;i++)
    {
        //cout<<"jet "<<i<<"\n";
        const GenJet& jet = (*m_vbfJets)[i];

        // jet requirements
        if(fabs(jet.Eta())<1.5) continue;
        if(fabs(jet.Eta())>3.) continue;
        if(jet.Pt()<30.) continue;

        vector<const Tower*> l1jets         = matchJet(jet);
        vector<const Tower*> l1jetsHard     = matchJetHard(jet);
        vector<const Tower*> l1jetsHardNoTh = matchJetHardNoThreshold(jet);

        const Tower* l1jet         = (l1jets.size()>0 ? l1jets[0] : 0);
        const Tower* l1jetHard     = (l1jetsHard.size()>0 ? l1jetsHard[0] : 0);
        const Tower* l1jetHardNoTh = (l1jetsHardNoTh.size()>0 ? l1jetsHardNoTh[0] : 0);

        if(!l1jet || !l1jetHard || !l1jetHardNoTh) continue;
        if(l1jet->calibratedEt()==0. || l1jetHard->calibratedEt()==0. || l1jetHardNoTh->calibratedEt()==0.) continue;
        // match with core
        const Tower* l1jetCore         = matchJetCore(*l1jet);
        const Tower* l1jetCoreHard     = matchJetCoreHard(*l1jetHard);
        const Tower* l1jetCoreHardNoTh = matchJetCoreHardNoThreshold(*l1jetHardNoTh);
        if(!l1jetCore || !l1jetCoreHard || !l1jetCoreHardNoTh) continue;
        if(l1jetCore->calibratedEt()==0. || l1jetCoreHard->calibratedEt()==0. || l1jetCoreHardNoTh->calibratedEt()==0.) continue;

        // match with 0.4 
        const Tower* l1jet04         = matchJet04(*l1jet);
        const Tower* l1jet04Hard     = matchJet04Hard(*l1jetHard);
        const Tower* l1jet04HardNoTh = matchJet04HardNoThreshold(*l1jetHardNoTh);
        if(!l1jet04 || !l1jet04Hard || !l1jet04HardNoTh) continue;
        if(l1jet04->calibratedEt()==0. || l1jet04Hard->calibratedEt()==0. || l1jet04HardNoTh->calibratedEt()==0.) continue;

        //cout<<"  ET = "<<l1jet->calibratedEt()<<"\n";

        if(l1jet->calibratedEt()<5.) continue;

        m_histos.FillHisto(10+hoffset, jet.Et(), weight, sysNum);
        m_histos.FillHisto(11+hoffset, l1jet->calibratedEt(), weight, sysNum);
        m_histos.FillHisto(12+hoffset, jet.E(), weight, sysNum);
        m_histos.FillHisto(13+hoffset, l1jet->calibratedEnergy(), weight, sysNum);

        // checks
        double responseThPU     = (l1jet->calibratedEt() - l1jetHardNoTh->calibratedEt())/l1jetHardNoTh->calibratedEt();
        double responseTrueThPU = (l1jet->calibratedEt() - jet.Pt())/jet.Pt();
        double responseTrueTh   = (l1jetHard->calibratedEt() - jet.Pt())/jet.Pt();
        double responseTrue     = (l1jetHardNoTh->calibratedEt() - jet.Pt())/jet.Pt();

        double coreFraction     = l1jetCore->calibratedEt()/l1jet->calibratedEt();
        double coreDeta         = l1jetCore->eta() - l1jet->eta();
        double coreDphi         = TVector2::Phi_mpi_pi(l1jetCore->phi() - l1jet->phi());
        double coreDR           = sqrt(coreDeta*coreDeta + coreDphi*coreDphi);

        double jet04Fraction     = l1jet04->calibratedEt()/l1jet->calibratedEt();
        double jet04Deta         = l1jet04->eta() - l1jet->eta();
        double jet04Dphi         = TVector2::Phi_mpi_pi(l1jet04->phi() - l1jet->phi());
        double jet04DR           = sqrt(jet04Deta*jet04Deta + jet04Dphi*jet04Dphi);

        m_histos.FillHisto(50+hoffset, responseThPU    , weight, sysNum);
        m_histos.FillHisto(51+hoffset, responseTrueThPU, weight, sysNum);
        m_histos.FillHisto(52+hoffset, responseTrueTh  , weight, sysNum);
        m_histos.FillHisto(53+hoffset, responseTrue    , weight, sysNum);

        m_histos.FillHisto(60+hoffset, coreFraction, weight, sysNum);
        m_histos.FillHisto(61+hoffset, coreDR, weight, sysNum);

        m_histos.FillHisto(70+hoffset, jet04Fraction, weight, sysNum);
        m_histos.FillHisto(71+hoffset, jet04DR, weight, sysNum);


        m_histos.Fill1BinHisto(100+hoffset, jet.Et(), l1jet->calibratedEt()/jet.Et(), weight, sysNum);
        m_histos.Fill1BinHisto(120+hoffset, jet.Et(), l1jet->calibratedEt(), weight, sysNum);
        m_histos.Fill1BinHisto(140+hoffset, jet.E(), l1jet->calibratedEnergy()/jet.E(), weight, sysNum);
        m_histos.Fill1BinHisto(160+hoffset, jet.E(), l1jet->calibratedEnergy(), weight, sysNum);

        /////////////////////////////////////////////////
        /// Regression tree filling
        /////////////////////////////////////////////////
        // regression targets
        double diffPU           = (l1jetHard->calibratedEnergy() - l1jet->calibratedEnergy())/(double)l1jet->nHits();
        double diffPUCore       = (l1jetCoreHard->calibratedEnergy() - l1jetCore->calibratedEnergy())/(double)l1jetCore->nHits();
        double diffPU04         = (l1jet04Hard->calibratedEnergy() - l1jet04->calibratedEnergy())/(double)l1jet04->nHits();
        //double correctionThHard = l1jetHardNoTh->calibratedEt()/l1jetHard->calibratedEt();
        double correctionTh     = l1jetHardNoTh->calibratedEt()/l1jet->calibratedEt();
        double correctionThCore = l1jetCoreHardNoTh->calibratedEt()/l1jetCore->calibratedEt();
        double correctionTh04   = l1jet04HardNoTh->calibratedEt()/l1jet04->calibratedEt();
        double correctionTruth  = jet.Et()/l1jet->calibratedEt();
        double correctionInvTruth  = l1jet->calibratedEt()/jet.Et();
        double correctionETruth  = jet.E()/l1jet->calibratedEnergy();
        double correctionEInvTruth  = l1jet->calibratedEnergy()/jet.E();

        // jet energy / core energy
        double jetOverCore     = l1jet->calibratedEt()/l1jetCore->calibratedEt();
        double jetOver04       = l1jet->calibratedEt()/l1jet04->calibratedEt();

        // compute NHits in 5 first ECAL layers in trigger region
        int triggerRegion = m_jetAlgo.triggerRegionIndex(l1jet->eta(), l1jet->hits()[0]->zside(), l1jet->hits()[0]->sector(), l1jet->hits()[0]->subsector());
        int nhits = 0;
        nhits += m_jetAlgo.triggerRegionHits(triggerRegion, 1, 3);
        nhits += m_jetAlgo.triggerRegionHits(triggerRegion, 2, 3);
        nhits += m_jetAlgo.triggerRegionHits(triggerRegion, 3, 3);
        nhits += m_jetAlgo.triggerRegionHits(triggerRegion, 4, 3);
        nhits += m_jetAlgo.triggerRegionHits(triggerRegion, 5, 3);

        // fill calibration tree
        m_histos.FillHisto(400+hoffset, diffPU, weight, sysNum);
        m_histos.FillHisto(401+hoffset, diffPUCore, weight, sysNum);
        m_histos.FillHisto(402+hoffset, diffPU04, weight, sysNum);
        m_histos.FillHisto(403+hoffset, correctionTh, weight, sysNum);
        m_histos.FillHisto(404+hoffset, correctionThCore, weight, sysNum);
        m_histos.FillHisto(405+hoffset, correctionTh04, weight, sysNum);
        m_histos.FillHisto(406+hoffset, correctionTruth, weight, sysNum);
        m_histos.FillHisto(407+hoffset, correctionInvTruth, weight, sysNum);
        m_histos.FillHisto(408+hoffset, correctionETruth, weight, sysNum);
        m_histos.FillHisto(409+hoffset, correctionEInvTruth, weight, sysNum);
        m_histos.FillHisto(410+hoffset, l1jet->eta(), weight, sysNum);
        m_histos.FillHisto(411+hoffset, l1jet->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(412+hoffset, jet.E(), weight, sysNum);
        m_histos.FillHisto(413+hoffset, l1jet->calibratedEt(), weight, sysNum);
        m_histos.FillHisto(414+hoffset, jet.Et(), weight, sysNum);
        m_histos.FillHisto(415+hoffset, nhits, weight, sysNum);
        m_histos.FillHisto(416+hoffset, jetOverCore, weight, sysNum);
        m_histos.FillHisto(417+hoffset, jetOver04, weight, sysNum);

        m_histos.FillNtuple(500+hoffset, event().run(), event().event(), weight, sysNum);

    }




}



/*****************************************************************/
void AnalysisJetCalibration::fillHistosTaus()
/*****************************************************************/
{

    //short sysNum = 0;
    //float weight = 1.;
    //int hoffset  = 3000;


}




/*****************************************************************/
void AnalysisJetCalibration::selectMaxMassJetPair()
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
void AnalysisJetCalibration::selectVBFJetPair()
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
bool AnalysisJetCalibration::overlapWithLeptonOrNeutrino(const GenJet& jet)
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
vector<const Tower*> AnalysisJetCalibration::matchJet(const TLorentzVector& jet)
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
vector<const Tower*> AnalysisJetCalibration::matchJetHard(const TLorentzVector& jet)
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
        if(dr<0.2) jets.push_back(&cl);
    }
    sort(jets.begin(),jets.end(), AnHiMa::towerSort);
    return jets;
}

/*****************************************************************/
vector<const Tower*> AnalysisJetCalibration::matchJetHardNoThreshold(const TLorentzVector& jet)
/*****************************************************************/
{
    vector<const Tower*> jets;
    double eta = jet.Eta();
    double phi = jet.Phi();
    // with PU
    for(const auto& cl : m_jetsHardNoThreshold)
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
const Tower* AnalysisJetCalibration::matchJetCore(const Tower& l1jet)
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

/*****************************************************************/
const Tower* AnalysisJetCalibration::matchJetCoreHard(const Tower& l1jet)
/*****************************************************************/
{
    const Tower* core = 0;
    double drMin = 99999.;
    double eta = l1jet.eta();
    double phi = l1jet.phi();
    for(const auto& c : m_jetCoresHard)
    {
        double deta = (c.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(c.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.05 && dr<drMin)
        {
            if(core) cout<<"WARNING: L1 hard jet matched to several jet cores\n";
            core = &c;
            drMin = dr;
        }
    }
    if(!core)
    {
        cout<<"WARNING: cannot find hard core jet\n";
    }
    return core;
}

/*****************************************************************/
const Tower* AnalysisJetCalibration::matchJetCoreHardNoThreshold(const Tower& l1jet)
/*****************************************************************/
{
    const Tower* core = 0;
    double drMin = 99999.;
    double eta = l1jet.eta();
    double phi = l1jet.phi();
    for(const auto& c : m_jetCoresHardNoThreshold)
    {
        double deta = (c.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(c.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.05 && dr<drMin)
        {
            if(core) cout<<"WARNING: L1 hard no th jet matched to several jet cores\n";
            core = &c;
            drMin = dr;
        }
    }
    if(!core)
    {
        cout<<"WARNING: cannot find hard no th core jet\n";
    }
    return core;
}

/*****************************************************************/
const Tower* AnalysisJetCalibration::matchJet04(const Tower& l1jet)
/*****************************************************************/
{
    const Tower* jet04 = 0;
    double drMin = 99999.;
    double eta = l1jet.eta();
    double phi = l1jet.phi();
    for(const auto& c : m_jets04)
    {
        double deta = (c.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(c.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.1 && dr<drMin)
        {
            if(jet04) cout<<"WARNING: L1 jet matched to several 0.4 jets\n";
            jet04 = &c;
            drMin = dr;
        }
    }
    if(!jet04)
    {
        cout<<"WARNING: cannot find 0.4 jet\n";
    }
    return jet04;
}

/*****************************************************************/
const Tower* AnalysisJetCalibration::matchJet04Hard(const Tower& l1jet)
/*****************************************************************/
{
    const Tower* jet04 = 0;
    double drMin = 99999.;
    double eta = l1jet.eta();
    double phi = l1jet.phi();
    for(const auto& c : m_jets04Hard)
    {
        double deta = (c.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(c.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.1 && dr<drMin)
        {
            if(jet04) cout<<"WARNING: L1 hard jet matched to several 0.4 jets\n";
            jet04 = &c;
            drMin = dr;
        }
    }
    if(!jet04)
    {
        cout<<"WARNING: cannot find hard 0.4 jet\n";
    }
    return jet04;
}

/*****************************************************************/
const Tower* AnalysisJetCalibration::matchJet04HardNoThreshold(const Tower& l1jet)
/*****************************************************************/
{
    const Tower* jet04 = 0;
    double drMin = 99999.;
    double eta = l1jet.eta();
    double phi = l1jet.phi();
    for(const auto& c : m_jets04HardNoThreshold)
    {
        double deta = (c.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(c.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.1 && dr<drMin)
        {
            if(jet04) cout<<"WARNING: L1 hard no th jet matched to several 0.4 jets\n";
            jet04 = &c;
            drMin = dr;
        }
    }
    if(!jet04)
    {
        cout<<"WARNING: cannot find hard no th 0.4 jet\n";
    }
    return jet04;
}

