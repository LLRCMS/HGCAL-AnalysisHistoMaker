/**
 *  @file  SimHit.h
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



#ifndef SIMHIT_H
#define SIMHIT_H

namespace AnHiMa
{
  class SimHit
  {
    public:
      SimHit();
      SimHit(int,int,int,int,float,float,float,float);
      SimHit(int,int,int,int,int,float,float,float,float,float,float);
      ~SimHit();


      int detid() const {return m_detid;}
      int type() const {return m_type;}
      int layer() const {return m_layer;}
      int sector() const {return m_sector;}
      int bin() const {return m_bin;}
      float x() const {return m_x;}
      float y() const {return m_y;}
      float z() const {return m_z;}
      float eta() const {return m_eta;}
      float phi() const {return m_phi;}
      float energy() const {return m_energy;}

    private:
      int   m_detid;
      int   m_type;
      int   m_layer;
      int   m_sector;
      int   m_bin;
      float m_x;
      float m_y;
      float m_z;
      float m_eta;
      float m_phi;
      float m_energy;

  };
};


#endif
