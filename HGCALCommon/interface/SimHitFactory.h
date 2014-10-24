



#ifndef SIMHITFACTORY_H
#define SIMHITFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

namespace AnHiMa
{
    typedef std::vector<SimHit> SimHitCollection;


    class SimHitFactory : public IObjectFactory<SimHitCollection>
    {
        public:
            SimHitFactory() {};
            ~SimHitFactory() {};

            virtual void initialize(IEvent*, TChain*);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);

            int                    m_hit_n;
            std::vector<unsigned> *m_hit_detid;
            std::vector<int>      *m_hit_subdet;
            std::vector<int>      *m_hit_cell;
            std::vector<int>      *m_hit_sector;
            std::vector<int>      *m_hit_subsector;
            std::vector<int>      *m_hit_layer;
            std::vector<int>      *m_hit_zside;
            std::vector<float>    *m_hit_energy;
            std::vector<float>    *m_hit_eta;
            std::vector<float>    *m_hit_phi;
            std::vector<float>    *m_hit_x;
            std::vector<float>    *m_hit_y;
            std::vector<float>    *m_hit_z;


    };
}


#endif
