/**
 *  @file  GenJet.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    07/22/2014
 *
 *  @internal
 *     Created :  07/22/2014
 * Last update :  07/22/2014 04:38:58 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef GENJET_H
#define GENJET_H

#include "TLorentzVector.h"

namespace AnHiMa
{
    class GenJet: public TLorentzVector
    {
        public:
            GenJet();
            ~GenJet();

            void setEmEnergy (float energy)        {m_emEnergy     = energy;}
            void setHadEnergy (float energy)       {m_hadEnergy     = energy;}
            void setInvisibleEnergy (float energy) {m_invisibleEnergy     = energy;}
            void setPuIndex (int ipu) {m_ipu     = ipu;}


            float emEnergy()        const {return m_emEnergy;}
            float hadEnergy()       const {return m_hadEnergy;}
            float invisibleEnergy() const {return m_invisibleEnergy;}
            float puIndex()         const {return m_ipu;}


        private:
            float m_emEnergy;
            float m_hadEnergy;
            float m_invisibleEnergy;
            float m_ipu;
    };
    bool genJetSort(const GenJet* jet1, const GenJet* jet2);

};


#endif
