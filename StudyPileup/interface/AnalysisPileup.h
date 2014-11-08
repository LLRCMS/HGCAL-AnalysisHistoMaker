/**
 *  @file  AnalysisPileup.h
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




#ifndef ANALYSISPILEUP_H
#define ANALYSISPILEUP_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TH2F.h>

#include <string>
#include <vector>




namespace AnHiMa
{

    class AnalysisPileup: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisPileup();
            ~AnalysisPileup();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            

    };
}


#endif
