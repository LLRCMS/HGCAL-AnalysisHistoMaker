/**
 *  @file  ProjectionMapProducer.h
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


#ifndef PROJECTIONMAPPRODUCER
#define PROJECTIONMAPPRODUCER

#include <vector>
#include <map>
#include <algorithm>

#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TGraph.h"

#include "AnHiMaHGCAL/HGCALGeometry/interface/HGCALNavigator.h"


namespace AnHiMa
{

    class ProjectionMapProducer
    {
        public:
            ProjectionMapProducer();
            ~ProjectionMapProducer();
            void initialize();
            void produce();

        private:
            void findBinning();

            TFile* m_output;
            TH2I* m_projectedBins;
            TH2I* m_unProjected;

            std::map< std::pair<unsigned,unsigned>, unsigned> m_cellToBin;
            std::vector< std::vector< std::pair<unsigned,unsigned> > > m_binToCell;

            int m_xN, m_yN;
            double m_xMin, m_xMax, m_yMin, m_yMax;


            HGCALNavigator m_hgcalNavigator;


            static const int NCELLS;

    };

};


#endif
