/**
 *  @file  AnalysisIdentification.h
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
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TH2F.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>




namespace AnHiMa
{

    class AnalysisIdentification: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisIdentification();
            ~AnalysisIdentification();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            void matchSuperCluster();

            TriggerEGammaAlgorithm m_egammaAlgo;


            std::vector<Tower> m_seeds;
            std::vector<Tower> m_clusters;
            std::vector<SuperCluster> m_superClusters;

            const GenParticle* m_genParticle;
            const SuperCluster* m_superCluster;

            std::string m_sample;
            float m_etCut;

    };
}


#endif
