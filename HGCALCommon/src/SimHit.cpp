

#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

#include "TMath.h"

#include <math.h>
#include <iostream>

using namespace AnHiMa;
using namespace std;

/*****************************************************************/
SimHit::SimHit():
  m_detid(0),
  m_type(0),
  m_layer(0),
  m_sector(0),
  m_bin(0),
  m_x(0.),
  m_y(0.),
  m_z(0.),
  m_eta(0.),
  m_phi(0.),
  m_energy(0.)
/*****************************************************************/
{
}

/*****************************************************************/
SimHit::SimHit(int type, int layer, int sector, int bin,
    float x, float y, float z, float energy):
  m_type(type),
  m_layer(layer),
  m_sector(sector),
  m_bin(bin),
  m_x(x),
  m_y(y),
  m_z(z),
  m_energy(energy)
/*****************************************************************/
{
  // compute eta and phi
  double rho = sqrt((double)x*(double)x+(double)y*(double)y+(double)z*(double)z);
  double phi = TMath::ATan2((double)y,(double)x);
  double eta = 0;
  if (rho>=z) eta=0.5*TMath::Log( (rho+(double)z)/(rho-(double)z) );
  m_eta = eta;
  m_phi = phi;

  // compute detector id
  int detid = 0;
  detid |= bin; // 13 bits
  detid |= (sector<<13); // 6 bits
  int abslayer = (layer>=0 ? layer : layer+63);
  detid |= (abslayer<<19); // 5 bits
  detid |= (type<<24); // 2 bits
  m_detid = detid;

}


/*****************************************************************/
SimHit::SimHit(int detid, int type, int layer, int sector, int bin,
    float x, float y, float z, 
    float eta, float phi, float energy):
  m_detid(detid),
  m_type(type),
  m_layer(layer),
  m_sector(sector),
  m_bin(bin),
  m_x(x),
  m_y(y),
  m_z(z),
  m_eta(eta),
  m_phi(phi),
  m_energy(energy)
/*****************************************************************/
{
}



/*****************************************************************/
SimHit::~SimHit()
/*****************************************************************/
{
}
