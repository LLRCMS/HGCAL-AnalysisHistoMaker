/**
 *  @file  AnalysisThreshold.h
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




#ifndef ANALYSISTHRESHOLD_H
#define ANALYSISTHRESHOLD_H

#include "AnHiMaHGCAL/Core/interface/IAnalysis.h"
#include "AnHiMaHGCAL/Core/interface/EventAware.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"

#include <TRandom3.h>
#include <TLorentzVector.h>
#include <TH2F.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>




namespace AnHiMa
{

    class AnalysisThreshold: public IAnalysis, EventAware<EventHGCAL>
    {
        public:
            AnalysisThreshold();
            ~AnalysisThreshold();

            bool initialize(const std::string& parameterFile);

            void execute();

        private:
            void fillHistos();
            void fillPileupParams();
            void fillRegionHits();
            void findMaxHit(int i);
            void findElectronHits();
            void findElectronHitsBelowThreshold();
            float pileupThreshold(float eta, int layer, int nhits);
            int triggerRegionIndex(float eta, int zside, int sector, int subsector);
            int triggerRegionHits(int triggerRegion, int layer);

            const SimHit* m_chosenHit;
            const GenParticle* m_genParticle;

            std::string m_pileupParamsFile;
            std::map< std::pair<int,int>, std::pair<float,float> > m_pileupParams;
            std::map< std::pair<int,int>, std::vector< std::pair<HGCEEDetId, float> > > m_regionSimHits; // detid (with eta) in trigger regions and layers inside a trigger region
            std::map< std::pair<int,int>, std::vector< const SimHit* > > m_regionHitsAboveTh;
            double m_electronHitEnergy;
            std::vector<HGCEEDetId> m_electronHits;
            std::vector<const SimHit*> m_electronHitsBelowThreshold;


            

    };
}


#endif
