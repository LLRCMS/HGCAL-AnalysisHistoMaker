



#ifndef PROJECTEDSIMHITFACTORY_H
#define PROJECTEDSIMHITFACTORY_H

#include <vector>
#include <string>
#include <map>

#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"


namespace AnHiMa
{
    class EventHGCAL;

    typedef std::vector< SimHit > SimHitSingleCollection;

    class ProjectedSimHitFactory : public IObjectFactory<SimHitSingleCollection>
    {
        public:
            enum HitType
            {
                ALL     = 0,
                HARDINT = 1
            };

        public:
            ProjectedSimHitFactory() {};
            ~ProjectedSimHitFactory() {};

            virtual void initialize(IEvent*, TChain*);
            void initialize(IEvent*, TChain*, const EventHGCAL* eventHGCAL, const std::string& projectionMappingFile, HitType);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);
            void fillProjectionMapping(const std::string& mappingFile);
            HGCEEDetId projectCell(const HGCEEDetId& cellid);

            HitType m_type;
            const EventHGCAL* m_event;
            TowerCalibrator m_calibrator;

            // x-y projection mapping
            std::map< std::pair<unsigned,unsigned>, unsigned> m_projectionMapping;
    };
}


#endif
