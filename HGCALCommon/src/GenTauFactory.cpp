/**
 *  @file  GenTauFactory.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    11/01/2015
 *
 *  @internal
 *     Created :  11/01/2015
 * Last update :  11/01/2015 16:07:25
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */







#include "AnHiMaHGCAL/HGCALCommon/interface/GenTauFactory.h"

using namespace AnHiMa;
using namespace std;



/*****************************************************************/
void GenTauFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, GenTauFactory::callback);
    connectVariables(chain);

}



/*****************************************************************/
void GenTauFactory::connectVariables(TChain* chain)
/*****************************************************************/
{
    m_gentau_n         = 0;
    m_gentau_energy    = 0;
    m_gentau_pt        = 0;
    m_gentau_eta       = 0;
    m_gentau_phi       = 0;
    m_gentau_decay     = 0;


    chain->SetBranchStatus("gentau_n"         , true);
    chain->SetBranchStatus("gentau_energy"    , true);
    chain->SetBranchStatus("gentau_pt"        , true);
    chain->SetBranchStatus("gentau_eta"       , true);
    chain->SetBranchStatus("gentau_phi"       , true);
    chain->SetBranchStatus("gentau_decay"     , true);


    chain->SetBranchAddress("gentau_n"         , &m_gentau_n);
    chain->SetBranchAddress("gentau_energy"    , &m_gentau_energy );
    chain->SetBranchAddress("gentau_pt"        , &m_gentau_pt     );
    chain->SetBranchAddress("gentau_eta"       , &m_gentau_eta    );
    chain->SetBranchAddress("gentau_phi"       , &m_gentau_phi    );
    chain->SetBranchAddress("gentau_decay"     , &m_gentau_decay  );

}



/*****************************************************************/
void GenTauFactory::callback(void* object)
/*****************************************************************/
{
    GenTauFactory* myself = reinterpret_cast<GenTauFactory*>(object);
    myself->update();
}


/*****************************************************************/
void GenTauFactory::update()
/*****************************************************************/
{
    m_data.clear();
    for(int i=0;i<m_gentau_n;i++)
    {
        GenTau gen;
        gen.setDecay     ( (*m_gentau_decay)[i] );
        gen.SetPtEtaPhiE ( (*m_gentau_pt)[i], (*m_gentau_eta)[i], (*m_gentau_phi)[i], (*m_gentau_energy)[i] );

        m_data.push_back( gen );
    }
}
