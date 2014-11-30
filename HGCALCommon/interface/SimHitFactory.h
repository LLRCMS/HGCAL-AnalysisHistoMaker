



#ifndef SIMHITFACTORY_H
#define SIMHITFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

namespace AnHiMa
{
    typedef std::vector<SimHit> SimHitCollection;

    bool simHitSort(const SimHit* hit1, const SimHit* hit2);


    class SimHitFactory : public IObjectFactory<SimHitCollection>
    {
        public:
            enum HitType
            {
                ALL     = 0,
                HARDINT = 1
            };

        public:
            SimHitFactory() {};
            ~SimHitFactory() {};

            virtual void initialize(IEvent*, TChain*);
            void initialize(IEvent*, TChain*, HitType);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);

            HitType m_type;

            int                    m_hit_n;
            std::vector<uint32_t> *m_hit_detid;
            std::vector<uint8_t>  *m_hit_subdet;
            std::vector<uint16_t> *m_hit_cell;
            std::vector<uint8_t>  *m_hit_sector;
            std::vector<char>     *m_hit_subsector;
            std::vector<uint8_t>  *m_hit_layer;
            std::vector<char>     *m_hit_zside;
            std::vector<float>    *m_hit_energy;
            std::vector<float>    *m_hit_eta;
            std::vector<float>    *m_hit_phi;
            std::vector<float>    *m_hit_x;
            std::vector<float>    *m_hit_y;
            std::vector<float>    *m_hit_z;


    };
}


#endif
