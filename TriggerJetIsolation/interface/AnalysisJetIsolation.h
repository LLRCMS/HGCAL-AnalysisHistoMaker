/**
 *  @file  AnalysisJetIsolation.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    07/02/2015
 *
 *  @internal
 *     Created :  07/02/2015
 * Last update :  07/02/2015 14:20:42
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */






#ifndef ANALYSISJETISOLATION_H
#define ANALYSISJETISOLATION_H

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

    class AnalysisJetIsolation: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisJetIsolation();
            ~AnalysisJetIsolation();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            // histo filling
            void fillHistos();
            void fillHistosVBFJets(int hoffset=0);
            void fillHistosMinBiasJets(int hoffset=0);

            // Utility functions
            void selectMaxMassJetPair();
            void selectVBFJetPair();
            bool overlapWithLeptonOrNeutrino(const GenJet& jet);
            std::pair<int, const Tower*> matchJet(const TLorentzVector& jet, const std::vector<Tower>& l1jets);
            double truncatedArea(double etaJ, double phiJ, double R, double eta0);


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
            std::map<std::string, std::vector<Tower> > m_jets;

            // jet selection
            TParticlePair<GenJet>* m_maxMassJets;
            TParticlePair<GenJet>* m_vbfJets;

            TRandom3 m_random;

            std::string m_jetType;

    };
}


#endif
