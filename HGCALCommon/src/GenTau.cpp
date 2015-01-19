/**
 *  @file  GenTau.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/01/2015
 *
 *  @internal
 *     Created :  11/01/2015
 * Last update :  11/01/2015 16:07:32
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include "AnHiMaHGCAL/HGCALCommon/interface/GenTau.h"

#include "TMath.h"

#include <math.h>
#include <iostream>

using namespace AnHiMa;
using namespace std;

/*****************************************************************/
GenTau::GenTau():TLorentzVector(),
  m_decay(0)
/*****************************************************************/
{
}




/*****************************************************************/
GenTau::~GenTau()
/*****************************************************************/
{
}
