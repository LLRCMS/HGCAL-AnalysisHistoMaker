/**
 *  @file  GenJetFactory.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    22/10/2014
 *
 *  @internal
 *     Created :  22/10/2014
 * Last update :  22/10/2014 21:31:42
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#include "AnHiMaHGCAL/HGCALCommon/interface/GenJetFactory.h"

using namespace AnHiMa;
using namespace std;



/*****************************************************************/
void GenJetFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, GenJetFactory::callback);
    connectVariables(chain);

}



/*****************************************************************/
void GenJetFactory::connectVariables(TChain* chain)
/*****************************************************************/
{
    m_genjet_n         = 0;
    m_genjet_energy    = 0;
    m_genjet_emenergy  = 0;
    m_genjet_hadenergy = 0;
    m_genjet_invenergy = 0;
    m_genjet_pt        = 0;
    m_genjet_eta       = 0;
    m_genjet_phi       = 0;


    chain->SetBranchStatus("genjet_n"         , true);
    chain->SetBranchStatus("genjet_energy"    , true);
    chain->SetBranchStatus("genjet_emenergy"  , true);
    chain->SetBranchStatus("genjet_hadenergy" , true);
    chain->SetBranchStatus("genjet_invenergy" , true);
    chain->SetBranchStatus("genjet_pt"        , true);
    chain->SetBranchStatus("genjet_eta"       , true);
    chain->SetBranchStatus("genjet_phi"       , true);


    chain->SetBranchAddress("genjet_n"         , &m_genjet_n);
    chain->SetBranchAddress("genjet_energy"    , &m_genjet_energy );
    chain->SetBranchAddress("genjet_emenergy"  , &m_genjet_emenergy );
    chain->SetBranchAddress("genjet_hadenergy" , &m_genjet_hadenergy );
    chain->SetBranchAddress("genjet_invenergy" , &m_genjet_invenergy );
    chain->SetBranchAddress("genjet_pt"        , &m_genjet_pt     );
    chain->SetBranchAddress("genjet_eta"       , &m_genjet_eta    );
    chain->SetBranchAddress("genjet_phi"       , &m_genjet_phi    );

}



/*****************************************************************/
void GenJetFactory::callback(void* object)
/*****************************************************************/
{
    GenJetFactory* myself = reinterpret_cast<GenJetFactory*>(object);
    myself->update();
}


/*****************************************************************/
void GenJetFactory::update()
/*****************************************************************/
{
    m_data.clear();
    for(int i=0;i<m_genjet_n;i++)
    {
        GenJet gen;
        gen.setEmEnergy        ( (*m_genjet_emenergy)[i] );
        gen.setHadEnergy       ( (*m_genjet_hadenergy)[i] );
        gen.setInvisibleEnergy ( (*m_genjet_invenergy)[i] );
        gen.SetPtEtaPhiE ( (*m_genjet_pt)[i], (*m_genjet_eta)[i], (*m_genjet_phi)[i], (*m_genjet_energy)[i] );

        m_data.push_back( gen );
    }
}
