/**
 *  @file  IEvent.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/21/2012
 *
 *  @internal
 *     Created :  12/21/2012
 * Last update :  12/21/2012 02:56:20 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <vector>
#include <TChain.h>
#include <TBranch.h>

#include "AnHiMaHGCAL/Core/interface/IEvent.h"

using namespace AnHiMa;
using namespace std;

/*****************************************************************/
IEvent::IEvent():
    m_inputChain(NULL)
/*****************************************************************/
{
}



/*****************************************************************/
IEvent::~IEvent()
/*****************************************************************/
{
}

/*****************************************************************/
void IEvent::connectVariables(TChain* inputChain)
/*****************************************************************/
{
    m_inputChain = inputChain;
    inputChain->SetBranchStatus("*", false);
    
}

/*****************************************************************/
void IEvent::registerCallback(void* objectPtr, CallBack function)
/*****************************************************************/
{
    m_listOfCallBacks.push_back( std::make_pair(objectPtr, function) );
}


/*****************************************************************/
void IEvent::update()
/*****************************************************************/
{
    for(unsigned int fu = 0; fu<m_listOfCallBacks.size();fu++)
    {
        (m_listOfCallBacks[fu].second)(m_listOfCallBacks[fu].first);
    }
}
