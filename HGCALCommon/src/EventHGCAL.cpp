/**
 *  @file  EventHGCAL.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    05/29/2014
 *
 *  @internal
 *     Created :  05/29/2014
 * Last update :  05/29/2014 05:17:36 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#include <iostream>
#include "AnHiMaHGCAL/HGCALCommon/interface/EventHGCAL.h"
#include <TChain.h>
#include "TMath.h"


using namespace AnHiMa;
using namespace std;

/*****************************************************************/
EventHGCAL::EventHGCAL():IEvent()
/*****************************************************************/
{


}



/*****************************************************************/
EventHGCAL::~EventHGCAL()
/*****************************************************************/
{
}


/*****************************************************************/
void EventHGCAL::connectVariables(TChain* inputChain)
/*****************************************************************/
{
    IEvent::connectVariables(inputChain);

    // Branches to read
    inputChain->SetBranchStatus("run"       , true);
    inputChain->SetBranchStatus("lumi"      , true);
    inputChain->SetBranchStatus("event"     , true);
    inputChain->SetBranchStatus("npu"       , true);


    // Connect branches
    inputChain->SetBranchAddress("run"       , &m_run);
    inputChain->SetBranchAddress("lumi"      , &m_lumi);
    inputChain->SetBranchAddress("event"     , &m_event);
    inputChain->SetBranchAddress("npu"       , &m_npu);


    m_simhitFactory.initialize(this, inputChain);

}


/*****************************************************************/
bool EventHGCAL::passSelection(int selection)
/*****************************************************************/
{
    return true;
}


