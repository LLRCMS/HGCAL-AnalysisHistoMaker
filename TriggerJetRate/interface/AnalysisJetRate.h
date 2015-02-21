/**
 *  @file  AnalysisJetRate.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    21/02/2015
 *
 *  @internal
 *     Created :  21/02/2015
 * Last update :  21/02/2015 13:41:37
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */






#ifndef ANALYSISJETRATE_H
#define ANALYSISJETRATE_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
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

    class AnalysisJetRate: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisJetRate();
            ~AnalysisJetRate();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();

            TriggerJetAlgorithm m_jetAlgo;
            TriggerEGammaAlgorithm m_egammaAlgo;
            std::vector<Tower> m_emSeeds;
            std::vector<Tower> m_emClusters;
            std::vector<Tower> m_hadClusters;
            std::vector<Tower> m_clusters;
            std::vector<SuperCluster> m_superClusters;
            std::map<std::string, std::vector<Tower> > m_jets;
            std::vector<const Tower*> m_sortedJets;


    };
}


#endif
