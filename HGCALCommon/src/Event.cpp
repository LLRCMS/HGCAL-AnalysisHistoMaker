/**
 *  @file  Event.cpp
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
#include "AnHiMaHGCAL/HGCALCommon/interface/Event.h"
#include <TChain.h>
#include "TMath.h"


using namespace AnHiMa;
using namespace std;

/*****************************************************************/
Event::Event():IEvent(),
    m_energySum(0.)
/*****************************************************************/
{


}



/*****************************************************************/
Event::~Event()
/*****************************************************************/
{
}


/*****************************************************************/
void Event::connectVariables(TChain* inputChain)
/*****************************************************************/
{
    IEvent::connectVariables(inputChain);

    // Branches to read
    inputChain->SetBranchStatus("run"       , true);
    inputChain->SetBranchStatus("lumi"      , true);
    inputChain->SetBranchStatus("event"     , true);

    inputChain->SetBranchStatus("ngen"      , true);
    inputChain->SetBranchStatus("gen_id"    , true);
    inputChain->SetBranchStatus("gen_status", true);
    inputChain->SetBranchStatus("gen_pt"    , true);
    inputChain->SetBranchStatus("gen_eta"   , true);
    inputChain->SetBranchStatus("gen_phi"   , true);
    inputChain->SetBranchStatus("gen_en"    , true);


    inputChain->SetBranchStatus("nhits"     , true);
    inputChain->SetBranchStatus("hit_type"  , true);
    inputChain->SetBranchStatus("hit_layer" , true);
    inputChain->SetBranchStatus("hit_sec"   , true);
    inputChain->SetBranchStatus("hit_bin"   , true);
    inputChain->SetBranchStatus("hit_x"     , true);
    inputChain->SetBranchStatus("hit_y"     , true);
    inputChain->SetBranchStatus("hit_z"     , true);
    inputChain->SetBranchStatus("hit_edep"  , true);

    inputChain->SetBranchStatus("hard_nhits"     , true);
    inputChain->SetBranchStatus("hard_hit_type"  , true);
    inputChain->SetBranchStatus("hard_hit_layer" , true);
    inputChain->SetBranchStatus("hard_hit_sec"   , true);
    inputChain->SetBranchStatus("hard_hit_bin"   , true);
    inputChain->SetBranchStatus("hard_hit_x"     , true);
    inputChain->SetBranchStatus("hard_hit_y"     , true);
    inputChain->SetBranchStatus("hard_hit_z"     , true);
    inputChain->SetBranchStatus("hard_hit_edep"  , true);



    // Connect branches
    inputChain->SetBranchAddress("run"       , &run);
    inputChain->SetBranchAddress("lumi"      , &lumi);
    inputChain->SetBranchAddress("event"     , &event);

    inputChain->SetBranchAddress("ngen"      , &ngen);
    inputChain->SetBranchAddress("gen_id"    , gen_id);
    inputChain->SetBranchAddress("gen_status", gen_status);
    inputChain->SetBranchAddress("gen_pt"    , gen_pt);
    inputChain->SetBranchAddress("gen_eta"   , gen_eta);
    inputChain->SetBranchAddress("gen_phi"   , gen_phi);
    inputChain->SetBranchAddress("gen_en"    , gen_en);

    inputChain->SetBranchAddress("nhits"     , &nhits);
    inputChain->SetBranchAddress("hit_type"  , hit_type);
    inputChain->SetBranchAddress("hit_layer" , hit_layer);
    inputChain->SetBranchAddress("hit_sec"   , hit_sec);
    inputChain->SetBranchAddress("hit_bin"   , hit_bin);
    inputChain->SetBranchAddress("hit_x"     , hit_x);
    inputChain->SetBranchAddress("hit_y"     , hit_y);
    inputChain->SetBranchAddress("hit_z"     , hit_z);
    inputChain->SetBranchAddress("hit_edep"  , hit_edep);

    inputChain->SetBranchAddress("hard_nhits"     , &hard_nhits);
    inputChain->SetBranchAddress("hard_hit_type"  , hard_hit_type);
    inputChain->SetBranchAddress("hard_hit_layer" , hard_hit_layer);
    inputChain->SetBranchAddress("hard_hit_sec"   , hard_hit_sec);
    inputChain->SetBranchAddress("hard_hit_bin"   , hard_hit_bin);
    inputChain->SetBranchAddress("hard_hit_x"     , hard_hit_x);
    inputChain->SetBranchAddress("hard_hit_y"     , hard_hit_y);
    inputChain->SetBranchAddress("hard_hit_z"     , hard_hit_z);
    inputChain->SetBranchAddress("hard_hit_edep"  , hard_hit_edep);


}


/*****************************************************************/
bool Event::passSelection(int selection)
/*****************************************************************/
{
    return true;
}



/*****************************************************************/
void Event::buildObjects()
/*****************************************************************/
{
    m_energySum = 0.;
    m_simHits.clear();
    m_simHits.reserve(nhits);
    for(int i=0;i<nhits;i++)
    {
        // compute eta and phi
        //double x = hit_x[i];
        //double y = hit_y[i];
        //double z = hit_z[i];
        //double rho = sqrt((double)x*(double)x+(double)y*(double)y+(double)z*(double)z);
        //double phi = TMath::ATan2((double)y,(double)x);
        //double eta = 0;
        //if (rho>=z) eta=0.5*TMath::Log( (rho+(double)z)/(rho-(double)z) );
        //hit_eta[i] = eta;
        //hit_phi[i] = phi;


        if(hit_type[i]!=0) continue;
        if(hit_layer[i]==0) continue;
        m_simHits.push_back(SimHit(hit_type[i],hit_layer[i],hit_sec[i],hit_bin[i],hit_x[i],hit_y[i],hit_z[i],hit_edep[i]));
        m_energySum += hit_edep[i];
    }

    m_hardEnergySum = 0.;
    m_hardSimHits.clear();
    m_hardSimHits.reserve(hard_nhits);
    for(int i=0;i<hard_nhits;i++)
    {

        if(hard_hit_type[i]!=0) continue;
        if(hard_hit_layer[i]==0) continue;
        m_hardSimHits.push_back(SimHit(hard_hit_type[i],hard_hit_layer[i],hard_hit_sec[i],hard_hit_bin[i],hard_hit_x[i],hard_hit_y[i],hard_hit_z[i],hard_hit_edep[i]));
        m_hardEnergySum += hard_hit_edep[i];
    }
}

