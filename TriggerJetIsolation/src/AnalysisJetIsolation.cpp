/**
 *  @file  AnalysisJetIsolation.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    07/02/2015
 *
 *  @internal
 *     Created :  07/02/2015
 * Last update :  07/02/2015 14:28:32
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */








#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/TriggerJetIsolation/interface/AnalysisJetIsolation.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisJetIsolation::AnalysisJetIsolation():IAnalysis(),
    m_maxMassJets(0),
    m_vbfJets(0),
    m_random(561034),
    m_jetType("VBF")
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisJetIsolation::~AnalysisJetIsolation()
/*****************************************************************/
{
    if(m_maxMassJets) delete m_maxMassJets;
}



/*****************************************************************/
bool AnalysisJetIsolation::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_egammaAlgo.initialize(event(), m_reader.params());
    m_jetAlgo.initialize(event(), m_reader.params());

    m_jetType = m_reader.params().GetValue("JetType", "VBF");

    m_jets["Size01"] = vector<Tower>();
    m_jets["Size02"] = vector<Tower>();
    m_jets["Size04"] = vector<Tower>();
    m_jets["Size05"] = vector<Tower>();

    cout<<"INFO: JetType = "<<m_jetType<<"\n";

    return true;
}


/*****************************************************************/
void AnalysisJetIsolation::execute()
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
    for(auto& col_jets : m_jets)
    {
        col_jets.second.clear();
    }

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
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets["Size01"], 0.1); 
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets["Size02"], 0.2); 
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets["Size04"], 0.4);
    m_jetAlgo.coneClustering(event(), simhits, m_superClusters, m_jets["Size05"], 0.5);


    if(m_jetType=="VBF")
    {
        selectMaxMassJetPair();
        selectVBFJetPair();
    }

    fillHistos();
}



/*****************************************************************/
void AnalysisJetIsolation::fillHistos()
/*****************************************************************/
{
    if(m_jetType=="VBF") fillHistosVBFJets(0);
    else if(m_jetType=="MinBias") fillHistosMinBiasJets(0);
    m_jetAlgo.pileupSubtraction(m_jets["Size01"], "0.1");
    m_jetAlgo.pileupSubtraction(m_jets["Size02"], "0.2");
    m_jetAlgo.pileupSubtraction(m_jets["Size04"], "0.4");
    m_jetAlgo.pileupSubtraction(m_jets["Size05"], "0.5");
    if(m_jetType=="VBF") fillHistosVBFJets(1000);
    else if(m_jetType=="MinBias") fillHistosMinBiasJets(1000);
    m_jetAlgo.thresholdCorrection(m_jets["Size01"], "0.1");
    m_jetAlgo.thresholdCorrection(m_jets["Size02"], "0.2");
    m_jetAlgo.thresholdCorrection(m_jets["Size04"], "0.4");
    m_jetAlgo.thresholdCorrection(m_jets["Size05"], "0.5");
    if(m_jetType=="VBF") fillHistosVBFJets(2000);
    else if(m_jetType=="MinBias") fillHistosMinBiasJets(2000);
}



/*****************************************************************/
void AnalysisJetIsolation::fillHistosVBFJets(int hoffset)
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;


    if(!m_vbfJets) return;
    for(unsigned i=0;i<2;i++)
    {
        const GenJet& jet = (*m_vbfJets)[i];

        // jet requirements
        if(fabs(jet.Eta())<1.5) continue;
        if(fabs(jet.Eta())>3.) continue;
        if(jet.Pt()<10.) continue;

        pair<unsigned, const Tower*> id_l1jet02 = matchJet(jet, m_jets["Size02"]);

        const Tower* l1jet02 = id_l1jet02.second;
        int jetindex = id_l1jet02.first;

        if(!l1jet02) continue;
        if(l1jet02->calibratedEt()<15.) continue;

        const SuperCluster* roi = &(m_superClusters[jetindex]);
        //
        const Tower* l1jet01 = &(m_jets["Size01"][jetindex]);
        const Tower* l1jet04 = &(m_jets["Size04"][jetindex]);
        const Tower* l1jet05 = &(m_jets["Size05"][jetindex]);


        // compute areas of cones
        unsigned nbins = 20;
        double binwidth = 0.05;
        vector<double> drBins(nbins,0.);
        vector<double> areas(nbins,0.);
        for(unsigned int drBin=0; drBin<nbins; drBin++)
        {
            drBins[drBin] = (double)drBin*binwidth+(binwidth/2.);
            double R = ((double)drBin+1.)*binwidth;
            double etaJ = fabs(roi->weightedEta())-2.25; // take the ROI center, and shift the endcaps around 0
            double phiJ = roi->weightedPhi();
            areas[drBin] = truncatedArea(etaJ, phiJ, R, 0.75);
        }
        vector<double> sumEBins = m_jetAlgo.computeJetProfile(*l1jet04, roi->weightedEta(), roi->weightedPhi());

        // check distance between jets of different sizes
        double deta01         = l1jet01->eta() - l1jet02->eta();
        double deta04         = l1jet04->eta() - l1jet02->eta();
        double deta05         = l1jet05->eta() - l1jet02->eta();
        double dphi01         = TVector2::Phi_mpi_pi(l1jet01->phi() - l1jet02->phi());
        double dphi04         = TVector2::Phi_mpi_pi(l1jet04->phi() - l1jet02->phi());
        double dphi05         = TVector2::Phi_mpi_pi(l1jet05->phi() - l1jet02->phi());
        double dR01           = sqrt(deta01*deta01 + dphi01*dphi01);
        double dR04           = sqrt(deta04*deta04 + dphi04*dphi04);
        double dR05           = sqrt(deta05*deta05 + dphi05*dphi05);
        //
        double area01 = areas[9];
        double area02 = areas[19];
        double area04 = areas[39];



        m_histos.FillHisto(10+hoffset, dR01, weight, sysNum);
        m_histos.FillHisto(11+hoffset, dR04, weight, sysNum);
        m_histos.FillHisto(12+hoffset, dR05, weight, sysNum);


        double previousArea = 0.;
        for(unsigned int drBin=0;drBin<nbins;drBin++)
        {
            double deltaR = drBins[drBin];
            double sumE = sumEBins[drBin];
            double area = areas[drBin] - previousArea;
            previousArea = areas[drBin];
            m_histos.Fill2BinHisto(50+hoffset, l1jet02->calibratedEt(), deltaR, sumE/area, weight, sysNum);
        }


        // fill isolation tree
        m_histos.FillHisto(500+hoffset,  roi->weightedEta(), weight, sysNum);
        m_histos.FillHisto(501+hoffset, l1jet02->eta(), weight, sysNum);
        m_histos.FillHisto(502+hoffset, l1jet02->calibratedEt(), weight, sysNum);
        m_histos.FillHisto(503+hoffset, l1jet02->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(504+hoffset, l1jet01->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(505+hoffset, l1jet04->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(506+hoffset, l1jet05->calibratedEnergy(), weight, sysNum);
        //
        m_histos.FillHisto(511+hoffset, area01, weight, sysNum);
        m_histos.FillHisto(512+hoffset, area02, weight, sysNum);
        m_histos.FillHisto(513+hoffset, area04, weight, sysNum);


        m_histos.FillNtuple(600+hoffset, event().run(), event().event(), weight, sysNum);

    }




}


/*****************************************************************/
void AnalysisJetIsolation::fillHistosMinBiasJets(int hoffset)
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;

    for(unsigned i=0; i<m_jets["Size02"].size(); i++)
    {
        const Tower* l1jet02 = &(m_jets["Size02"][i]);
        int jetindex = i;

        if(l1jet02->calibratedEt()<15.) continue;

        const SuperCluster* roi = &(m_superClusters[jetindex]);
        //
        const Tower* l1jet01 = &(m_jets["Size01"][jetindex]);
        const Tower* l1jet04 = &(m_jets["Size04"][jetindex]);
        const Tower* l1jet05 = &(m_jets["Size05"][jetindex]);

        // compute areas of cones
        unsigned nbins = 20;
        double binwidth = 0.05;
        vector<double> drBins(nbins,0.);
        vector<double> areas(nbins,0.);
        for(unsigned int drBin=0; drBin<nbins; drBin++)
        {
            drBins[drBin] = (double)drBin*binwidth+(binwidth/2.);
            double R = ((double)drBin+1.)*binwidth;
            double etaJ = fabs(roi->weightedEta())-2.25; // take the ROI center, and shift the endcaps around 0
            double phiJ = roi->weightedPhi();
            areas[drBin] = truncatedArea(etaJ, phiJ, R, 0.75);
        }
        vector<double> sumEBins = m_jetAlgo.computeJetProfile(*l1jet04, roi->weightedEta(), roi->weightedPhi());


        // check distance between jets of different sizes
        double deta01         = l1jet01->eta() - l1jet02->eta();
        double deta04         = l1jet04->eta() - l1jet02->eta();
        double deta05         = l1jet05->eta() - l1jet02->eta();
        double dphi01         = TVector2::Phi_mpi_pi(l1jet01->phi() - l1jet02->phi());
        double dphi04         = TVector2::Phi_mpi_pi(l1jet04->phi() - l1jet02->phi());
        double dphi05         = TVector2::Phi_mpi_pi(l1jet05->phi() - l1jet02->phi());
        double dR01           = sqrt(deta01*deta01 + dphi01*dphi01);
        double dR04           = sqrt(deta04*deta04 + dphi04*dphi04);
        double dR05           = sqrt(deta05*deta05 + dphi05*dphi05);
        //
        double area01 = areas[9];
        double area02 = areas[19];
        double area04 = areas[39];


        m_histos.FillHisto(10+hoffset, dR01, weight, sysNum);
        m_histos.FillHisto(11+hoffset, dR04, weight, sysNum);
        m_histos.FillHisto(12+hoffset, dR05, weight, sysNum);

        double previousArea = 0.;
        for(unsigned int drBin=0;drBin<nbins;drBin++)
        {
            double deltaR = drBins[drBin];
            double sumE = sumEBins[drBin];
            double area = areas[drBin] - previousArea;
            previousArea = areas[drBin];
            m_histos.Fill2BinHisto(50+hoffset, l1jet02->calibratedEt(), deltaR, sumE/area, weight, sysNum);
        }


        // fill isolation tree
        m_histos.FillHisto(500+hoffset,  roi->weightedEta(), weight, sysNum);
        m_histos.FillHisto(501+hoffset, l1jet02->eta(), weight, sysNum);
        m_histos.FillHisto(502+hoffset, l1jet02->calibratedEt(), weight, sysNum);
        m_histos.FillHisto(503+hoffset, l1jet02->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(504+hoffset, l1jet01->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(505+hoffset, l1jet04->calibratedEnergy(), weight, sysNum);
        m_histos.FillHisto(506+hoffset, l1jet05->calibratedEnergy(), weight, sysNum);
        //
        m_histos.FillHisto(511+hoffset, area01, weight, sysNum);
        m_histos.FillHisto(512+hoffset, area02, weight, sysNum);
        m_histos.FillHisto(513+hoffset, area04, weight, sysNum);

        m_histos.FillNtuple(600+hoffset, event().run(), event().event(), weight, sysNum);

    }




}




/*****************************************************************/
void AnalysisJetIsolation::selectMaxMassJetPair()
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
void AnalysisJetIsolation::selectVBFJetPair()
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
bool AnalysisJetIsolation::overlapWithLeptonOrNeutrino(const GenJet& jet)
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
pair<int, const Tower*> AnalysisJetIsolation::matchJet(const TLorentzVector& jet, const vector<Tower>& l1jets)
/*****************************************************************/
{
    pair<int, const Tower*> l1jet = make_pair<int, const Tower*>(-1, 0);
    double eta = jet.Eta();
    double phi = jet.Phi();
    double drMin = 9999.;
    // with PU
    for(unsigned i=0; i<l1jets.size(); i++)
    {
        const auto& j = l1jets[i];
        double deta = (j.eta() - eta);
        double dphi = TVector2::Phi_mpi_pi(j.phi() - phi);
        double dr = sqrt(deta*deta + dphi*dphi);
        if(dr<0.2 && dr<drMin)
        {
            l1jet = make_pair((int)i,&j);
            drMin = dr;
        }
    }
    return l1jet;
}


/*****************************************************************/
double AnalysisJetIsolation::truncatedArea(double etaJ, double phiJ, double R, double eta0)
/*****************************************************************/
{
    if(R>M_PI) cout<<"WARNING: AnalysisJetIsolation::truncatedArea(): R > Pi : circle will overlap itself\n";

    double deltaEta1 = fabs(etaJ-eta0);
    double deltaEta2 = fabs(etaJ+eta0);
    double deltaEta1Square = deltaEta1*deltaEta1;
    double deltaEta2Square = deltaEta2*deltaEta2;
    bool intersect1 = (R*R-deltaEta1Square>0 ? true : false);
    bool intersect2 = (R*R-deltaEta2Square>0 ? true : false);
    double phi11 = 0.;
    double phi12 = 0.;
    double phi21 = 0.;
    double phi22 = 0.;
    if((R<=deltaEta1 && intersect1) || (R<=deltaEta2 && intersect2))
        throw string("ERROR: AnalysisJetIsolation::truncatedArea(): Found intersections while it shouldn't");
    if((R>deltaEta1 && !intersect1) || (R>deltaEta2 && !intersect2))
        throw string("ERROR: AnalysisJetIsolation::truncatedArea(): Didn't found intersections while it should");

    if(intersect1)
    {
        phi11 = phiJ + sqrt(R*R-deltaEta1Square);
        phi12 = phiJ - sqrt(R*R-deltaEta1Square);
    }
    if(intersect2)
    {
        phi21 = phiJ + sqrt(R*R-deltaEta2Square);
        phi22 = phiJ - sqrt(R*R-deltaEta2Square);
    }

    double area = 0.;
    if(!intersect1 && !intersect2)
        area = M_PI*R*R;
    else if(intersect1 && !intersect2)
    {
        double deltaPhi = fabs(phi11-phi12);
        double alpha = 2.*asin(deltaPhi/(2.*R));
        double areaCircle = M_PI*R*R*(1.-alpha/(2.*M_PI));
        double areaTriangle = deltaEta1*deltaPhi/2.;
        area = areaCircle + areaTriangle;
    }
    else if(!intersect1 && intersect2)
    {
        double deltaPhi = fabs(phi21-phi22);
        double alpha = 2.*asin(deltaPhi/(2.*R));
        double areaCircle = M_PI*R*R*(1.-alpha/(2.*M_PI));
        double areaTriangle = deltaEta2*deltaPhi/2.;
        area = areaCircle + areaTriangle;
    }
    else if(intersect1 && intersect2)
    {
        double deltaPhi1 = fabs(phi11-phi12);
        double deltaPhi2 = fabs(phi21-phi22);
        double alpha1 = 2.*asin(deltaPhi1/(2.*R));
        double alpha2 = 2.*asin(deltaPhi2/(2.*R));
        double areaCircle = M_PI*R*R*(1.-(alpha1+alpha2)/(2.*M_PI));
        double areaTriangle = (deltaEta1*deltaPhi1 + deltaEta2*deltaPhi2)/2.;
        area = areaCircle + areaTriangle;
    }
    return area;

}
