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

namespace AnHiMa
{

    class TriggerEGammaAlgorithm
    {
        public:
            TriggerEGammaAlgorithm();
            ~TriggerEGammaAlgorithm();

            void initialize(const EventHGCAL& event, const std::string& pileupParamsFile);

            void seeding(const EventHGCAL& event, std::vector<Tower>& seeds);

        private:
            void fillPileupEstimators(const EventHGCAL& event);

            std::map< std::pair<int,int>, std::pair<float,float> > m_pileupParams;
            std::map< std::pair<int,int>, std::vector< std::pair<HGCEEDetId, float> > > m_regionSimHits; 
            std::map< std::pair<int,int>, std::vector< const SimHit* > > m_regionHitsAboveTh; 
                

    };

};



#endif
