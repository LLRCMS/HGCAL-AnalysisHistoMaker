



#ifndef GENJETFACTORY_H
#define GENJETFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/GenJet.h"

namespace AnHiMa
{
    typedef std::vector<GenJet> GenJetCollection;


    class GenJetFactory : public IObjectFactory<GenJetCollection>
    {
        public:
            GenJetFactory() {};
            ~GenJetFactory() {};

            virtual void initialize(IEvent*, TChain*);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);

            int                 m_genjet_n;
            std::vector<float> *m_genjet_energy;
            std::vector<float> *m_genjet_emenergy;
            std::vector<float> *m_genjet_hadenergy;
            std::vector<float> *m_genjet_invenergy;
            std::vector<float> *m_genjet_pt;
            std::vector<float> *m_genjet_eta;
            std::vector<float> *m_genjet_phi;
            std::vector<int>   *m_genjet_ipu;

    };
}


#endif
