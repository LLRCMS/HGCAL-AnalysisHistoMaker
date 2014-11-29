/**
 *  @file  TriggerEGammaAlgorithm.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    17/11/2014
 *
 *  @internal
 *     Created :  17/11/2014
 * Last update :  17/11/2014 18:51:26
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef TRIGGEREGAMMAALGORITHM_H
#define TRIGGEREGAMMAALGORITHM_H

#include <vector>
#include <string>

#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"
//
#include "AnHiMaHGCAL/HGCALCommon/interface/GBRForest.h"

namespace AnHiMa
{

    class TriggerEGammaAlgorithm
    {
        public:
            TriggerEGammaAlgorithm();
            ~TriggerEGammaAlgorithm();

            void initialize(const EventHGCAL& event, TEnv& params);

            void seeding(const EventHGCAL& event, std::vector<Tower>& seeds);
            void clustering(const EventHGCAL& event, const std::vector<Tower>& seeds, std::vector<Tower>& clusters, bool hardHits=false, bool pileupTh=true);
            void clusterPileupSubtraction(std::vector<Tower>& clusters);
            void clusterThresholdCorrection(std::vector<Tower>& clusters);
            void superClustering(const std::vector<Tower>& clusters, std::vector<SuperCluster>& superClusters);
            void superClusterCorrection(std::vector<SuperCluster>& superClusters);
            void idealClustering(const EventHGCAL& event, std::vector<Tower>& clusters, double eta0, double phi0);

            float pileupThreshold(float eta, int layer, int nhits);
            int triggerRegionIndex(float eta, int zside, int sector, int subsector);
            int triggerRegionHits(int triggerRegion, int layer);

        private:
            void fillPileupEstimators(const EventHGCAL& event);


            // pileup related 
            std::map< std::pair<int,int>, std::pair<float,float> > m_pileupParams;
            std::map< std::pair<int,int>, std::vector< std::pair<HGCEEDetId, float> > > m_regionSimHits; 
            std::map< std::pair<int,int>, std::vector< const SimHit* > > m_regionHitsAboveTh; 
            // clustering related
            std::map<int, float> m_clusterSizes; 
            // cluster energy corrections
            GBRForest* m_pileupSubtraction_up;
            GBRForest* m_pileupSubtraction_down;
            GBRForest* m_thresholdCorrection_up;
            GBRForest* m_thresholdCorrection_down;
            GBRForest* m_superClusterCorrection;
    };

};



#endif
