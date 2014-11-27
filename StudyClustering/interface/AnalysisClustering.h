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
            void fillHistos();
            void matchCluster(int i);
            void matchClusterHard(int i);
            void matchClusterHardNoPUThreshold(int i);
            void matchSuperCluster(int i);
            void matchSuperClusterHard(int i);
            void matchSuperClusterHardNoPUThreshold(int i);
            void matchSuperClusterCorr(int i);
            //
            const Tower* findSeedClusterHard(const SuperCluster* sc);
            const Tower* findSeedClusterHardNoPUThreshold(const SuperCluster* sc);

            TriggerEGammaAlgorithm m_egammaAlgo;


            std::vector<Tower> m_seeds;
            std::vector<Tower> m_clusters;
            std::vector<Tower> m_clustersHard;
            std::vector<Tower> m_clustersHardNoPUThreshold;
            std::vector<Tower> m_clustersCorr;
            std::vector<SuperCluster> m_superClusters;
            std::vector<SuperCluster> m_superClustersHard;
            std::vector<SuperCluster> m_superClustersHardNoPUThreshold;
            std::vector<SuperCluster> m_superClustersCorr;

            const GenParticle* m_genParticle;
            std::vector<const Tower*> m_matchedClusters;
            std::vector<const Tower*> m_matchedClustersHard;
            std::vector<const Tower*> m_matchedClustersHardNoPUThreshold;
            const SuperCluster* m_matchedSuperCluster;
            const SuperCluster* m_matchedSuperClusterHard;
            const SuperCluster* m_matchedSuperClusterHardNoPUThreshold;
            const SuperCluster* m_matchedSuperClusterCorr;


    };
}


#endif
