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
#include "AnHiMaHGCAL/HGCALCommon/interface/GenParticle.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHitFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TriggerHitFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/ProjectedSimHitFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TowerFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/GenParticleFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/GenJetFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/GenTauFactory.h"

#include "AnHiMaHGCAL/HGCALGeometry/interface/HGCALNavigator.h"

#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"

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
            void update();
            void initializeHitProjection(const std::string& projectionMappingFile);


            int event()  const {return m_event;}
            int lumi()   const {return m_lumi;}
            int run()    const {return m_run;}
            int npu()    const {return m_npu;}

            const HGCALNavigator& hgcalNavigator() const {return m_hgcalNavigator;}

            // hits
            const std::vector<SimHit>& simhits() const {return m_simhitFactoryAll.data()[0];}
            const std::vector<SimHit>& simhitsE() const {return m_simhitFactoryAll.data()[0];}
            const std::vector<SimHit>& simhitsFH() const {return m_simhitFactoryAll.data()[1];}
            const std::vector<SimHit>& simhitsBH() const {return m_simhitFactoryAll.data()[2];}
            const std::vector<SimHit>& hardsimhits() const {return m_simhitFactoryHard.data()[0];}
            const std::vector<SimHit>& hardsimhitsE() const {return m_simhitFactoryHard.data()[0];}
            const std::vector<SimHit>& hardsimhitsFH() const {return m_simhitFactoryHard.data()[1];}
            const std::vector<SimHit>& hardsimhitsBH() const {return m_simhitFactoryHard.data()[2];}
            const std::vector<const SimHit*>& sortedHits() const {return m_energySortedHits;}
            const std::vector<const SimHit*>& sortedSeedingHits() const {return m_energySortedSeedingHits;}
            //
            const std::vector<SimHit>& projectedhits() const {return m_projectedSimhitFactoryAll.data();}
            const std::vector<SimHit>& hardprojectedhits() const {return m_projectedSimhitFactoryHard.data();}
            //
            const std::vector<Tower>& towers() const {return m_towerFactory.data();}

            // gen info
            const std::vector<GenParticle>& genparticles() const {return m_genparticleFactory.data();}
            const std::vector<GenJet>& genjets() const {return m_genjetFactory.data();}
            const std::vector<GenTau>& gentaus() const {return m_gentauFactory.data();}

            const SimHit* simhit(int zside, int layer, int sector, int subsector, int cell) const;
            const SimHit* simhitFH(int zside, int layer, int sector, int subsector, int cell) const;
            const SimHit* simhitBH(int zside, int layer, int sector, int subsector, int cell) const;
            const SimHit* simhit(const HGCEEDetId& detid) const;
            const SimHit* simhitFH(const HGCHEDetId& detid) const;
            const SimHit* simhitBH(const HGCHEDetId& detid) const;
            const SimHit* hardsimhit(int zside, int layer, int sector, int subsector, int cell) const;
            const SimHit* hardsimhitFH(int zside, int layer, int sector, int subsector, int cell) const;
            const SimHit* hardsimhitBH(int zside, int layer, int sector, int subsector, int cell) const;
            const SimHit* hardsimhit(const HGCEEDetId& detid) const;
            const SimHit* hardsimhitFH(const HGCHEDetId& detid) const;
            const SimHit* hardsimhitBH(const HGCHEDetId& detid) const;

            unsigned cellIndex(int zside, int layer, int sector, int subsector, int cell) const;
            unsigned cellIndex(const HGCEEDetId& detid) const;
            unsigned cellIndex(const HGCHEDetId& detid) const;


        private:
            static void callback(void*);
            void buildHitsIndex();
            void buildHitsIndexFH();
            void buildHitsIndexBH();
            void sortHits();

            const double m_mip;


            int m_event;
            int m_lumi;
            int m_run;
            int m_npu;

            HGCALNavigator m_hgcalNavigator;

            SimHitFactory      m_simhitFactoryAll;
            SimHitFactory      m_simhitFactoryHard;
            TriggerHitFactory  m_triggerhitFactoryAll;
            TriggerHitFactory  m_triggerhitFactoryHard;
            ProjectedSimHitFactory m_projectedSimhitFactoryAll;
            ProjectedSimHitFactory m_projectedSimhitFactoryHard;
            TowerFactory       m_towerFactory;
            GenParticleFactory m_genparticleFactory;
            GenJetFactory      m_genjetFactory;
            GenTauFactory      m_gentauFactory;

            std::vector<const SimHit*> m_sortedHits;
            std::vector<const SimHit*> m_sortedHitsFH;
            std::vector<const SimHit*> m_sortedHitsBH;
            std::vector<const SimHit*> m_sortedHardHits;
            std::vector<const SimHit*> m_sortedHardHitsFH;
            std::vector<const SimHit*> m_sortedHardHitsBH;
            std::vector<const SimHit*> m_energySortedHits;
            std::vector<const SimHit*> m_energySortedSeedingHits;

    };
}

#endif
