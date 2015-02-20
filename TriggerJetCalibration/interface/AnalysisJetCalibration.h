/**
 *  @file  AnalysisJetCalibration.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    02/02/2015
 *
 *  @internal
 *     Created :  02/02/2015
 * Last update :  02/02/2015 11:22:13
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#ifndef ANALYSISJETCALIBRATION_H
#define ANALYSISJETCALIBRATION_H

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

    class AnalysisJetCalibration: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisJetCalibration();
            ~AnalysisJetCalibration();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            // histo filling
            void fillHistos();
            void fillHistosVBFJets(int hoffset=0);
            void fillHistosTaus();
            void fillHistosJet(const GenJet& jet, short sysNum, float weight, int hoffset);
            void fillHistosTau(const GenTau& tau, short sysNum, float weight, int hoffset);

            // Utility functions
            void selectMaxMassJetPair();
            void selectVBFJetPair();
            bool overlapWithLeptonOrNeutrino(const GenJet& jet);
            std::vector<const Tower*> matchJet(const TLorentzVector& jet);
            std::vector<const Tower*> matchJetHard(const TLorentzVector& jet);
            std::vector<const Tower*> matchJetHardNoThreshold(const TLorentzVector& jet);
            const Tower* matchJetCore(const Tower& l1jet);
            const Tower* matchJetCoreHard(const Tower& l1jet);
            const Tower* matchJetCoreHardNoThreshold(const Tower& l1jet);
            const Tower* matchJet04(const Tower& l1jet);
            const Tower* matchJet04Hard(const Tower& l1jet);
            const Tower* matchJet04HardNoThreshold(const Tower& l1jet);


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
            std::vector<Tower> m_jetCores;
            std::vector<Tower> m_jets04;
            std::vector<Tower> m_jetsHard;
            std::vector<Tower> m_jetCoresHard;
            std::vector<Tower> m_jets04Hard;
            std::vector<Tower> m_jetsHardNoThreshold;
            std::vector<Tower> m_jetCoresHardNoThreshold;
            std::vector<Tower> m_jets04HardNoThreshold;

            // jet selection
            TParticlePair<GenJet>* m_maxMassJets;
            TParticlePair<GenJet>* m_vbfJets;

            TRandom3 m_random;

    };
}


#endif
