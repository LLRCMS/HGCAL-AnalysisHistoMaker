



#ifndef GENPARTICLEFACTORY_H
#define GENPARTICLEFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/GenParticle.h"

namespace AnHiMa
{
    typedef std::vector<GenParticle> GenParticleCollection;


    class GenParticleFactory : public IObjectFactory<GenParticleCollection>
    {
        public:
            GenParticleFactory() {};
            ~GenParticleFactory() {};

            virtual void initialize(IEvent*, TChain*);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);

            int                m_gen_n;
            std::vector<int>   *m_gen_id;
            std::vector<int>   *m_gen_status;
            std::vector<float> *m_gen_energy;
            std::vector<float> *m_gen_pt;
            std::vector<float> *m_gen_eta;
            std::vector<float> *m_gen_phi;

    };
}


#endif
