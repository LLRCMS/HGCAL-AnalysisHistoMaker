



#ifndef TRIGGERHITFACTORY_H
#define TRIGGERHITFACTORY_H

#include <vector>
#include <array>

#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHitFactory.h"
#include "AnHiMaHGCAL/HGCALGeometry/interface/HGCALNavigator.h"

namespace AnHiMa
{
    struct SimplePosition
    {
        float x;
        float y;
    };

    class TriggerHitFactory : public IObjectFactory<SimHitCollection>
    {
        public:
            TriggerHitFactory() {};
            ~TriggerHitFactory() {};

            virtual void initialize(IEvent*, TChain*);
            void initialize(IEvent*, TChain*, const SimHitFactory& , const std::string&, const HGCALNavigator& );

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);
            unsigned cellIndex(int zside, int layer, int sector, int cell) const;

            std::array< std::map<short, short>, 30 > m_cellToTriggerCell;
            std::array< std::array< std::map<short, SimplePosition>, 30>, 18> m_triggerCellPositions;


            const SimHitFactory* m_hitFactory;

    };
}


#endif
