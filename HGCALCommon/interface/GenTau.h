/**
 *  @file  GenTau.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/01/2015
 *
 *  @internal
 *     Created :  11/01/2015
 * Last update :  11/01/2015 16:04:36
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#ifndef GENTAU_H
#define GENTAU_H

#include "TLorentzVector.h"

namespace AnHiMa
{
  class GenTau: public TLorentzVector
  {
    public:
      GenTau();
      ~GenTau();

      void setDecay (int decay)        {m_decay = decay;}


      int decay()        const {return m_decay;}


    private:
      float m_decay;
  };
};


#endif
