/**
 *  @file  GenParticle.h
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



#ifndef GENPARTICLE_H
#define GENPARTICLE_H

#include "TLorentzVector.h"

namespace AnHiMa
{
  class GenParticle: public TLorentzVector
  {
    public:
      GenParticle();
      ~GenParticle();

      void setId     (int id)      {m_id     = id;}
      void setStatus (int status)  {m_status    = status;}


      int  id()     const {return m_id;}
      int  status() const {return m_status;}


    private:
      int m_id;
      int m_status;

  };
};


#endif
