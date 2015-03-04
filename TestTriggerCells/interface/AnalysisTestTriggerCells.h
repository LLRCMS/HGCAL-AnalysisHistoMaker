/**
 *  @file  AnalysisTestTriggerCells.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    02/03/2015
 *
 *  @internal
 *     Created :  02/03/2015
 * Last update :  02/03/2015 11:21:34
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef ANALYSISTESTTRIGGERCELLS_H
#define ANALYSISTESTTRIGGERCELLS_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"

#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include <string>
#include <vector>
#include <map>
#include <algorithm>




namespace AnHiMa
{

    class AnalysisTestTriggerCells: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisTestTriggerCells();
            ~AnalysisTestTriggerCells();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            unsigned cellIndex(int zside, int layer, int sector, int cell) const;

            std::vector< std::map<short, short> > m_cellToTriggerCell;
            std::array< std::map<short, std::vector<short>>, 30> m_triggerCellToCell; 


    };
}


#endif
