/**
 *  @file  AnalysisProjection.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    23/01/2015
 *
 *  @internal
 *     Created :  23/01/2015
 * Last update :  23/01/2015 12:23:30
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */






#ifndef ANALYSISPROJECTION_H
#define ANALYSISPROJECTION_H

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
#include <TGraph.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>




namespace AnHiMa
{

    class AnalysisProjection: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisProjection();
            ~AnalysisProjection();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            void fillProjectionMapping(const std::string& inputFile);
            HGCEEDetId projectCell(const HGCEEDetId& cellid);

            int m_eventCount;
            std::map< std::pair<unsigned,unsigned>, unsigned> m_projectionMapping;

    };
}


#endif
