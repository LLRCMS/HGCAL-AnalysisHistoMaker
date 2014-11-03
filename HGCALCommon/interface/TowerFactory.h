



#ifndef TOWERFACTORY_H
#define TOWERFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"

namespace AnHiMa
{
    class EventHGCAL;

    typedef std::vector<Tower> TowerCollection;



    class TowerFactory : public IObjectFactory<TowerCollection>
    {
        public:
            enum HitType
            {
                ALL     = 0,
                HARDINT = 1
            };

        public:
            TowerFactory() {};
            ~TowerFactory() {};

            virtual void initialize(IEvent*, TChain*);
            void initialize(IEvent*, TChain*, EventHGCAL*);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);

            EventHGCAL* m_hgcalEvent;


    };
}


#endif
