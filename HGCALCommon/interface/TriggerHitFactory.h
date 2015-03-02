



#ifndef TRIGGERHITFACTORY_H
#define TRIGGERHITFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHitFactory.h"

namespace AnHiMa
{

    class TriggerHitFactory : public IObjectFactory<SimHitCollection>
    {
        public:
            TriggerHitFactory() {};
            ~TriggerHitFactory() {};

            virtual void initialize(IEvent*, TChain*);
            void initialize(IEvent*, TChain*, const SimHitFactory& , const std::string& );

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);
            unsigned cellIndex(int zside, int layer, int sector, int cell) const;

            std::vector< std::map<short, short> > m_cellToTriggerCell;

            const SimHitFactory* m_hitFactory;

    };
}


#endif
