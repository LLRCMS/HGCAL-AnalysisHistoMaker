/**
 *  @file  GenTauFactory.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/01/2015
 *
 *  @internal
 *     Created :  11/01/2015
 * Last update :  11/01/2015 16:06:01
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#ifndef GENTAUFACTORY_H
#define GENTAUFACTORY_H

#include <vector>
#include "AnHiMaHGCAL/Core/interface/IObjectFactory.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/GenTau.h"

namespace AnHiMa
{
    typedef std::vector<GenTau> GenTauCollection;


    class GenTauFactory : public IObjectFactory<GenTauCollection>
    {
        public:
            GenTauFactory() {};
            ~GenTauFactory() {};

            virtual void initialize(IEvent*, TChain*);

        private:
            virtual void connectVariables(TChain*);
            virtual void update();
            static void callback(void*);

            int                 m_gentau_n;
            std::vector<float> *m_gentau_energy;
            std::vector<float> *m_gentau_pt;
            std::vector<float> *m_gentau_eta;
            std::vector<float> *m_gentau_phi;
            std::vector<int>   *m_gentau_decay;

    };
}


#endif
