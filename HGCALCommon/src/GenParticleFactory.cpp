/**
 *  @file  GenParticleFactory.cpp
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





#include "AnHiMaHGCAL/HGCALCommon/interface/GenParticleFactory.h"

using namespace AnHiMa;
using namespace std;



/*****************************************************************/
void GenParticleFactory::initialize(IEvent* event, TChain* chain)
/*****************************************************************/
{
    event->registerCallback((void*)this, GenParticleFactory::callback);
    connectVariables(chain);

}



/*****************************************************************/
void GenParticleFactory::connectVariables(TChain* chain)
/*****************************************************************/
{
    m_gen_n      = 0;
    m_gen_id     = 0;
    m_gen_status = 0;
    m_gen_energy = 0;
    m_gen_pt     = 0;
    m_gen_eta    = 0;
    m_gen_phi    = 0;


    chain->SetBranchStatus("gen_n"      , true);
    chain->SetBranchStatus("gen_id"     , true);
    chain->SetBranchStatus("gen_status" , true);
    chain->SetBranchStatus("gen_energy" , true);
    chain->SetBranchStatus("gen_pt"     , true);
    chain->SetBranchStatus("gen_eta"    , true);
    chain->SetBranchStatus("gen_phi"    , true);


    chain->SetBranchAddress("gen_n"      , &m_gen_n);
    chain->SetBranchAddress("gen_id"     , &m_gen_id    );
    chain->SetBranchAddress("gen_status" , &m_gen_status);
    chain->SetBranchAddress("gen_energy" , &m_gen_energy );
    chain->SetBranchAddress("gen_pt"     , &m_gen_pt     );
    chain->SetBranchAddress("gen_eta"    , &m_gen_eta    );
    chain->SetBranchAddress("gen_phi"    , &m_gen_phi    );

}



/*****************************************************************/
void GenParticleFactory::callback(void* object)
/*****************************************************************/
{
    GenParticleFactory* myself = reinterpret_cast<GenParticleFactory*>(object);
    myself->update();
}


/*****************************************************************/
void GenParticleFactory::update()
/*****************************************************************/
{
    m_data.clear();
    for(int i=0;i<m_gen_n;i++)
    {
        GenParticle gen;
        gen.setId        ( (*m_gen_id)[i] );
        gen.setStatus    ( (*m_gen_status)[i] );
        gen.SetPtEtaPhiE ( (*m_gen_pt)[i], (*m_gen_eta)[i], (*m_gen_phi)[i], (*m_gen_energy)[i] );

        m_data.push_back( gen );
    }
}
