/**
 *  @file  NavigationTester.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/12/2014
 *
 *  @internal
 *     Created :  05/12/2014
 * Last update :  05/12/2014 16:57:43
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */


#ifndef NAVIGATIONTESTER
#define NAVIGATIONTESTER

#include <vector>
#include <map>
#include <algorithm>

#include "TFile.h"
#include "TGraph.h"

#include "AnHiMaHGCAL/HGCALGeometry/interface/HGCALNavigator.h"


namespace AnHiMa
{

    class NavigationTester
    {
        public:
            NavigationTester();
            ~NavigationTester();
            void initialize();
            void leftRightTest();

        private:
            HGCALNavigator m_hgcalNavigator;

            TFile* m_output;
            TGraph* m_cellsSource;
            TGraph* m_cells_1_0;
            TGraph* m_cells_2_0;
            TGraph* m_cells_3_0;
            TGraph* m_cells_m1_0;
            TGraph* m_cells_m2_0;
            TGraph* m_cells_m3_0;

    };

};


#endif
