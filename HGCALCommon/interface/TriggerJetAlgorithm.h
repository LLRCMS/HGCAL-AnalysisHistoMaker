/**
 *  @file  TriggerJetAlgorithm.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    14/01/2015
 *
 *  @internal
 *     Created :  14/01/2015
 * Last update :  14/01/2015 13:50:24
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#ifndef TRIGGERJETALGORITHM_H
#define TRIGGERJETALGORITHM_H

#include <vector>
#include <string>

#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SuperCluster.h"

namespace AnHiMa
{

    class TriggerJetAlgorithm
    {
        public:
            TriggerJetAlgorithm();
            ~TriggerJetAlgorithm();

            void initialize(const EventHGCAL& event, TEnv& params);

            void hcalSeeding(const EventHGCAL& event, const std::vector<SimHit>& hits, std::vector<Tower>& seeds, double clusterSize=3.);
            void superClustering(const EventHGCAL& event, const std::vector<Tower>& seeds, std::vector<SuperCluster>& superClusters, double coneSize=0.4);
            void coneClustering(const EventHGCAL& event, const std::vector<SimHit>& hits, const std::vector<SuperCluster>& superClusters, std::vector<Tower>& cones, double coneSize=0.4);

            float pileupThreshold(float eta, int layer, int nhits, int subdet=3);
            int triggerRegionIndex(float eta, int zside, int sector, int subsector);
            int triggerRegionHits(int triggerRegion, int layer, int subdet=3);

        private:
            void fillPileupEstimators(const EventHGCAL& event);

            // pileup related 
            std::map< std::pair<int,int>, std::pair<float,float> > m_pileupParamsECAL;
            std::map< std::pair<int,int>, std::pair<float,float> > m_pileupParamsHCAL;
            std::map< std::pair<int,int>, std::vector< std::pair<HGCEEDetId, float> > > m_regionSimHitsECAL; 
            std::map< std::pair<int,int>, std::vector< std::pair<HGCHEDetId, float> > > m_regionSimHitsHCAL;
            std::map< std::pair<int,int>, std::vector< const SimHit* > > m_regionHitsAboveThECAL; 
            std::map< std::pair<int,int>, std::vector< const SimHit* > > m_regionHitsAboveThHCAL; 
    };

};



#endif