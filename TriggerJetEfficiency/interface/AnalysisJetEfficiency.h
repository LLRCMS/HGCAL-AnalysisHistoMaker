/**
 *  @file  AnalysisJetEfficiency.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 13:40:35
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef ANALYSISJETEFFICIENCY_H
#define ANALYSISJETEFFICIENCY_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include "AnHiMaHGCAL/Core/interface/TParticlePair.h"

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

    class AnalysisJetEfficiency: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisJetEfficiency();
            ~AnalysisJetEfficiency();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();

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


            // L1 related
            TriggerEGammaAlgorithm m_egammaAlgo;
            TriggerJetAlgorithm m_jetAlgo;
            std::vector<Tower> m_emSeeds;
            std::vector<Tower> m_emClusters;
            std::vector<Tower> m_hadClusters;
            std::vector<Tower> m_clusters;
            std::vector<SuperCluster> m_superClusters;
            std::vector<Tower> m_jets;
            std::vector<Tower> m_jetCores;

            // jet selection
            TParticlePair<GenJet>* m_maxMassJets;
            TParticlePair<GenJet>* m_vbfJets;
            const GenJet* m_genJet;



            std::vector<int> m_thresholds;
            TTree** m_failAndPassTrees;
            double m_gen_pt[1];
            double m_gen_eta[1];
            int m_failAndPassBits[1][200];
            int m_failAndPassIDBits[1][200];


    };
}


#endif
