/**
 *  @file  AnalysisJets.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/01/2015
 *
 *  @internal
 *     Created :  05/01/2015
 * Last update :  05/01/2015 15:02:56
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef ANALYSISJETS_H
#define ANALYSISJETS_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/Core/interface/TParticlePair.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerEGammaAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerJetAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"


#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TH2F.h>
#include <TGraph.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>




namespace AnHiMa
{

    class AnalysisJets: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisJets();
            ~AnalysisJets();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            // histo filling
            void fillHistos();
            void fillHistosAllJets();
            void fillHistosVBFJets();
            void fillHistosRandom();
            void fillHistosTaus();
            void fillHistosJet(const GenJet& jet, short sysNum, float weight, int hoffset);
            void fillHistosTau(const GenTau& tau, short sysNum, float weight, int hoffset);

            // Utility functions
            void selectMaxMassJetPair();
            void selectVBFJetPair();
            bool overlapWithLeptonOrNeutrino(const GenJet& jet);
            std::vector<const Tower*> matchClusters(const TLorentzVector& jet);
            std::vector<const Tower*> matchClustersHard(const TLorentzVector& jet);
            std::vector<const Tower*> matchHadClusters(const TLorentzVector& jet);
            std::vector<const Tower*> matchHadClustersHard(const TLorentzVector& jet);
            std::vector<const SuperCluster*> matchROI(const TLorentzVector& jet);
            std::vector<const SuperCluster*> matchROIHard(const TLorentzVector& jet);
            std::vector<const Tower*> matchJet(const TLorentzVector& jet);
            std::vector<const Tower*> matchJetHard(const TLorentzVector& jet);
            Tower coneClustering(double eta0, double phi0, double coneSize=0.1, bool EH=false, bool hard=false);
            std::vector<Tower> coneClustering(double eta0, double phi0, const std::vector<double>& coneSizes, bool hard=false);

            // L1 related
            TriggerEGammaAlgorithm m_egammaAlgo;
            TriggerJetAlgorithm m_jetAlgo;
            std::vector<Tower> m_emSeeds;
            std::vector<Tower> m_emClusters;
            std::vector<Tower> m_emClustersHard;
            std::vector<Tower> m_hadClusters;
            std::vector<Tower> m_hadClustersHard;
            std::vector<Tower> m_clusters;
            std::vector<Tower> m_clustersHard;
            std::vector<SuperCluster> m_superClusters;
            std::vector<SuperCluster> m_superClustersHard;
            std::vector<Tower> m_jets;
            std::vector<Tower> m_jetsHard;

            std::map<int, TGraph*> m_bdtCuts;

            // jet selection
            TParticlePair<GenJet>* m_maxMassJets;
            TParticlePair<GenJet>* m_vbfJets;

            TRandom3 m_random;




    };
}


#endif
