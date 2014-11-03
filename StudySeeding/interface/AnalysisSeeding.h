/**
 *  @file  AnalysisSeeding.h
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




#ifndef ANALYSISSEEDING_H
#define ANALYSISSEEDING_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TH2F.h>

#include <string>
#include <vector>




namespace AnHiMa
{

    class AnalysisSeeding: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisSeeding();
            ~AnalysisSeeding();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos(int towerSize=1, bool cut=false);
            void buildTower(int size=1);
            std::vector<HGCEEDetId> choose2x2(const std::vector<HGCEEDetId>& detids);
            std::vector<HGCEEDetId> chooseNxN(const std::vector<HGCEEDetId>& detids, int N=2);
            void findMaxHit(int i);
            

            Tower m_tower;
            TowerCalibrator m_calibrator;
            const SimHit* m_chosenHit;
            const GenParticle* m_genParticle;

            TRandom3 m_random;

            std::string m_sample;

            std::vector<const SimHit*> m_usedHits;

    };
}


#endif
