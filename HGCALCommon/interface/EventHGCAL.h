/**
 *  @file  EventHGCAL.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/29/2014
 *
 *  @internal
 *     Created :  05/29/2014
 * Last update :  05/29/2014 04:45:27 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef EVENTHGCAL_H
#define EVENTHGCAL_H

#include "AnHiMaHGCAL/Core/interface/IEvent.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHitFactory.h"

#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <string>




namespace AnHiMa
{

    class EventHGCAL: public IEvent
    {
        public:
            EventHGCAL();
            ~EventHGCAL();

            bool passSelection(int selection=0);
            void connectVariables(TChain* inputChain);
            void buildObjects();

            int event()  const {return m_event;}
            int lumi()   const {return m_lumi;}
            int run()    const {return m_run;}
            int npu()    const {return m_npu;}

            const std::vector<SimHit>& simhits() const {return m_simhitFactory.data();}


        private:
            int m_event;
            int m_lumi;
            int m_run;
            int m_npu;

            SimHitFactory m_simhitFactory;

    };
}

#endif
