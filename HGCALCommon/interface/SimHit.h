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
      ~SimHit();

      void setDetid    (unsigned detid){m_detid     = detid;}
      void setSubdet   (int subdet)    {m_subdet    = subdet;}
      void setCell     (int cell)      {m_cell      = cell;}
      void setSector   (int sector)    {m_sector    = sector;}
      void setSubsector(int subsector) {m_subsector = subsector;}
      void setLayer    (int layer)     {m_layer     = layer;}
      void setZside    (int zside)     {m_zside     = zside;}
      void setEnergy   (float energy)  {m_energy    = energy;}
      void setEta      (float eta)     {m_eta       = eta;}
      void setPhi      (float phi)     {m_phi       = phi;}
      void setX        (float x)       {m_x         = x;}
      void setY        (float y)       {m_y         = y;}
      void setZ        (float z)       {m_z         = z;}


      unsigned detid()     const {return m_detid;}
      int      subdet()    const {return m_subdet;}
      int      cell()      const {return m_cell;}
      int      sector()    const {return m_sector;}
      int      subsector() const {return m_subsector;}
      int      layer()     const {return m_layer;}
      int      zside()     const {return m_zside;}
      float    energy()    const {return m_energy;}
      float    eta()       const {return m_eta;}
      float    phi()       const {return m_phi;}
      float    x()         const {return m_x;}
      float    y()         const {return m_y;}
      float    z()         const {return m_z;}


    private:
      unsigned m_detid;
      int      m_subdet;
      int      m_cell;
      int      m_sector;
      int      m_subsector;
      int      m_layer;
      int      m_zside;
      float    m_energy;
      float    m_eta;
      float    m_phi;
      float    m_x;
      float    m_y;
      float    m_z;

  };
};


#endif
