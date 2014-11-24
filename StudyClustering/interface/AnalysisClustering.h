/**
 *  @file  AnalysisClustering.h
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




#ifndef ANALYSISCLUSTERING_H
#define ANALYSISCLUSTERING_H

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

    class AnalysisClustering: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisClustering();
            ~AnalysisClustering();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos(bool fiducialCut);
            void matchCluster(int i);
            void matchSuperCluster(int i);

            TriggerEGammaAlgorithm m_egammaAlgo;


            std::vector<Tower> m_seeds;
            std::vector<Tower> m_clusters;
            std::vector<SuperCluster> m_superClusters;

            const GenParticle* m_genParticle;
            std::vector<const Tower*> m_matchedClusters;
            const SuperCluster* m_matchedSuperCluster;


    };
}


#endif
