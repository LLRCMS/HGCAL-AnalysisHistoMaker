

#include "AnHiMaHGCAL/HGCALCommon/interface/GenJet.h"

#include "TMath.h"

#include <math.h>
#include <iostream>

using namespace AnHiMa;
using namespace std;

/*****************************************************************/
GenJet::GenJet():TLorentzVector(),
  m_emEnergy(0),
  m_hadEnergy(0),
  m_invisibleEnergy(0)
/*****************************************************************/
{
}




/*****************************************************************/
GenJet::~GenJet()
/*****************************************************************/
{
}
