

#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

#include "TMath.h"

#include <math.h>
#include <iostream>

using namespace AnHiMa;
using namespace std;

/*****************************************************************/
SimHit::SimHit():
  m_detid(0),
  m_subdet(0),
  m_cell(0),
  m_sector(0),
  m_subsector(0),
  m_layer(0),
  m_zside(0),
  m_energy(0.),
  m_eta(0.),
  m_phi(0.),
  m_x(0.),
  m_y(0.),
  m_z(0.)
/*****************************************************************/
{
}




/*****************************************************************/
SimHit::~SimHit()
/*****************************************************************/
{
}
