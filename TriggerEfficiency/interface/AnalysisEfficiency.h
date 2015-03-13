/**
 *  @file  AnalysisEfficiency.h
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




#ifndef ANALYSISEFFICIENCY_H
#define ANALYSISEFFICIENCY_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerEGammaAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerJetAlgorithm.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

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

    class AnalysisEfficiency: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisEfficiency();
            ~AnalysisEfficiency();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            void matchSeed();
            void matchCluster();
            void matchSuperCluster();

            TriggerEGammaAlgorithm m_egammaAlgo;
            TriggerJetAlgorithm m_jetAlgo;


            std::vector<Tower> m_seeds;
            std::vector<Tower> m_clusters;
            std::vector<SuperCluster> m_superClusters;
            std::vector<Tower> m_idealClusters;
            std::vector<Tower> m_hcalCones;

            const GenParticle* m_genParticle;
            const Tower* m_matchedSeed;
            const Tower* m_matchedCluster;
            const SuperCluster* m_matchedSuperCluster;
            const Tower* m_matchedIdealCluster;
            const Tower* m_matchedHcalCone;

            std::vector<int> m_thresholds;
            TTree** m_failAndPassTrees;
            double m_gen_pt[1];
            double m_gen_eta[1];
            int m_failAndPassBits[1][100];
            int m_failAndPassIDBits[1][100];
            int m_failAndPassIDandHoEBits[1][100];

            std::map<int, TGraph*> m_bdtCuts;
            std::map<int, TGraph*> m_hoeCuts;


    };
}


#endif
