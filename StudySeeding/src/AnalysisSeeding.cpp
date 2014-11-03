/**
 *  @file  AnalysisSeeding.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    19/10/2014
 *
 *  @internal
 *     Created :  19/10/2014
 * Last update :  19/10/2014 13:44:32
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <iostream>

#include "TVector2.h"

#include "AnHiMaHGCAL/StudySeeding/interface/AnalysisSeeding.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"






using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisSeeding::AnalysisSeeding():IAnalysis(),
    m_chosenHit(0),
    m_genParticle(0),
    m_random(12345),
    m_sample("")
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisSeeding::~AnalysisSeeding()
/*****************************************************************/
{
}



/*****************************************************************/
bool AnalysisSeeding::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_sample = m_reader.params().GetValue("Sample", "minbias");



    return true;
}


/*****************************************************************/
void AnalysisSeeding::execute()
/*****************************************************************/
{
    event().update();
    if(!event().passSelection()) return;

    m_usedHits.clear();

    //unsigned n = (m_sample=="minbias" ? event().towers().size(): 2);
    unsigned n = (m_sample=="minbias" ? 500 : 2);
    //cout<<"Number of towers = "<<event().towers().size()<<"\n";

    vector<int> towerSizes = {1,2,4};
    for(unsigned i=0;i<n;i++)
    {
        findMaxHit(i);
        if(!m_chosenHit)
        {
            if(m_sample!="minbias") cout<<"[WARNING] Cannot match a hit close to the chosen (eta,phi) direction\n";
            continue;
        }
        for(auto isize=towerSizes.begin(); isize!=towerSizes.end(); isize++)
        {
            int size = *isize;
            buildTower(size);
            if(m_tower.energy()==0) 
            {
                cout<<"[WARNING] Tower has 0 energy\n";
                continue;
            }

            fillHistos(size, false);
            if(m_tower.calibratedEt()>2.) fillHistos(size, true);
        }
    }

}

/*****************************************************************/
void AnalysisSeeding::fillHistos(int towerSize, bool cut)
/*****************************************************************/
{
    //double mip = 0.000055;

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;
    switch(towerSize)
    {
        case 1: 
            hoffset = 0;
            break;
        case 2: 
            hoffset = 20000;
            break;
        case 4: 
            hoffset = 40000;
            break;
        default:
            break;
    }
    if(cut) hoffset += 10000;


    m_histos.FillHisto(0+hoffset, 0.5, weight, sysNum); // Number of events

    for(unsigned l=1; l<=30; l++)
    {
        double fraction = m_tower.layerEnergy(l)/m_tower.energy();
        m_histos.Fill2BinHisto(100+hoffset, fabs(m_tower.eta()), l, fraction, weight, sysNum);
        m_histos.Fill2BinHisto(200+hoffset, fabs(m_tower.eta()), l, m_tower.layerEnergy(l), weight, sysNum);
    }
    double energyFront  = m_tower.layersEnergy(1,10);
    double energyMiddle = m_tower.layersEnergy(11,20);
    double energyBack   = m_tower.layersEnergy(21,30);
    double fractionFront  = energyFront/m_tower.energy();
    double fractionMiddle = energyMiddle/m_tower.energy();
    double fractionBack   = energyBack/m_tower.energy();
    m_histos.Fill2BinHisto(300+hoffset, fabs(m_tower.eta()), 1, fractionFront, weight, sysNum);
    m_histos.Fill2BinHisto(300+hoffset, fabs(m_tower.eta()), 2, fractionMiddle, weight, sysNum);
    m_histos.Fill2BinHisto(300+hoffset, fabs(m_tower.eta()), 3, fractionBack, weight, sysNum);
    //
    m_histos.Fill2BinHisto(350+hoffset, fabs(m_tower.eta()), 1, energyFront, weight, sysNum);
    m_histos.Fill2BinHisto(350+hoffset, fabs(m_tower.eta()), 2, energyMiddle, weight, sysNum);
    m_histos.Fill2BinHisto(350+hoffset, fabs(m_tower.eta()), 3, energyBack, weight, sysNum);

    double middleOverFront = ( energyFront>0. ? energyMiddle/energyFront : 19.999);
    if(middleOverFront>=20.) middleOverFront = 19.999;
    m_histos.Fill1BinHisto(400+hoffset, fabs(m_tower.eta()), middleOverFront, weight, sysNum);

    m_histos.Fill1BinHisto(405+hoffset, fabs(m_tower.eta()), m_tower.energy(), weight, sysNum);
    m_histos.Fill1BinHisto(410+hoffset, fabs(m_tower.eta()), m_tower.et(), weight, sysNum);

    m_histos.Fill1BinHisto(420+hoffset, fabs(m_tower.eta()), m_tower.nHits(), weight, sysNum);
    m_histos.Fill1BinHisto(425+hoffset, fabs(m_tower.eta()), m_tower.nLayers(), weight, sysNum);
    m_histos.Fill1BinHisto(430+hoffset, fabs(m_tower.eta()), m_tower.longitudinalBarycenter(), weight, sysNum);

    m_histos.Fill1BinHisto(435+hoffset, fabs(m_tower.eta()), m_tower.calibratedEnergy(), weight, sysNum);
    m_histos.Fill1BinHisto(440+hoffset, fabs(m_tower.eta()), m_tower.calibratedEt(), weight, sysNum);

    // filling tree
    m_histos.FillHisto(1000+hoffset, m_tower.eta(), weight, sysNum);
    m_histos.FillHisto(1001+hoffset, m_tower.phi(), weight, sysNum);
    m_histos.FillHisto(1002+hoffset, m_tower.energy(), weight, sysNum);
    m_histos.FillHisto(1003+hoffset, m_tower.calibratedEnergy(), weight, sysNum);
    for(int l=1;l<=30;l++)
    {
        double fraction = m_tower.layerEnergy(l)/m_tower.energy();
        m_histos.FillHisto(1003+l+hoffset, fraction, weight, sysNum);
    }
    m_histos.FillHisto(1034+hoffset, m_tower.nHits(), weight, sysNum);
    m_histos.FillNtuple(1100+hoffset, event().run(), event().event(), weight, sysNum);
}


/*****************************************************************/
void AnalysisSeeding::buildTower(int size)
/*****************************************************************/
{
    double mip = 0.000055;
    m_tower.clear();

    HGCEEDetId detid(ForwardSubdetector(m_chosenHit->subdet()), m_chosenHit->zside(), m_chosenHit->layer(), m_chosenHit->sector(),m_chosenHit->subsector(), m_chosenHit->cell());
    double refEta = m_chosenHit->eta();
    double refPhi = m_chosenHit->phi();

    // find cells in the tower at (eta,phi)
    int layer = detid.layer();
    vector<HGCEEDetId> idLayers(31);
    idLayers[layer] = detid;
    //if(layer!=15)
    //{
    //    vector<HGCEEDetId> idLayer15 = event().hgcalNavigator().upProj(detid, 15-layer);
    //    if(idLayer15.size()>0) idLayers[15] = idLayer15[0];
    //    else
    //    {
    //        cout<<"[WARNING] Cannot find cell in layer 15\n";
    //        return;
    //    }
    //}
    for(int l=layer-1;l>0;l--)
    {
        vector<HGCEEDetId> ids;
        ids = event().hgcalNavigator().downProj(idLayers[l+1], 1, refEta, refPhi);
        if(ids.size()==0)
        {
            cout<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
            idLayers[l] = (HGCEEDetId)DetId(0);
        }
        else
        {
            idLayers[l] = ids[0];
        }
    }
    for(int l=layer+1;l<=30;l++)
    {
        vector<HGCEEDetId> ids;
        ids = event().hgcalNavigator().upProj(idLayers[l-1], 1, refEta, refPhi);
        if(ids.size()==0)
        {
            cout<<"[WARNING] Cannot find cell in layer "<<l<<"\n";
            idLayers[l] = (HGCEEDetId)DetId(0);
        }
        else
        {
            idLayers[l] = ids[0];
        }
    }

    //cout<<"Building 2x2 tower at "<<m_chosenHit->eta()<<", "<<m_chosenHit->phi()<<"\n";
    vector<HGCEEDetId> tower2x2 = chooseNxN(idLayers, size);
    //vector<HGCEEDetId> tower2x2_Bis = choose2x2(idLayers);

    // build tower from the list of cells
    //for(unsigned l=0; l<idLayers.size(); l++)
    for(auto itr=tower2x2.begin(); itr!=tower2x2.end(); itr++)
    {
        //HGCEEDetId id = idLayers[l];
        HGCEEDetId id = *itr;
        if(id==DetId(0)) continue;
        const SimHit* hit = event().simhit(id);
        //cout<<"  Possible hit in the tower: layer="<<id.layer()<<",cell="<<id.cell()<<"\n";
        if(!hit || hit->energy()<=mip) continue;
        m_tower.addHit(*hit);
    }
    m_calibrator.calibrate(m_tower);
    //cout<<" Tower energy  = "<<m_tower.energy()<<"\n";
    //cout<<" Tower NHits   = "<<m_tower.nHits()<<"\n";
    //cout<<" Tower NLayers = "<<m_tower.nLayers()<<"\n";



    //cout<<"Tower (eta,phi)=("<<m_tower.eta()<<","<<m_tower.phi()<<"): #hits="<<m_tower.nHits()<<", energy="<<m_tower.energy()<<"\n";



}

/*****************************************************************/
vector<HGCEEDetId> AnalysisSeeding::chooseNxN( const vector<HGCEEDetId>& detids, int N)
/*****************************************************************/
{
    vector<double> energy(4);
    energy[0] = 0.;
    energy[1] = 0.;
    energy[2] = 0.;
    energy[3] = 0.;

    vector< vector<DetId> > idsNeighbour(4);
    idsNeighbour[0] = vector<DetId>();
    idsNeighbour[1] = vector<DetId>();
    idsNeighbour[2] = vector<DetId>();
    idsNeighbour[3] = vector<DetId>();



    int n = N/2;
    int r = N%2;

    for(auto itr=detids.cbegin(); itr!=detids.end(); itr++)
    {
        const HGCEEDetId& detid = *itr;
        if(detid==DetId(0)) continue;

        // 
        for(int shiftx=-n; shiftx<=n; shiftx++)
        {
            for(int shifty=-n; shifty<=n; shifty++)
            {
                if(shiftx==0 && shifty==0) continue;

                DetId id = event().hgcalNavigator().topology()->offsetBy(detid, shiftx,  shifty);

                if(abs(shiftx)<n+r && abs(shifty)<n+r)
                {
                    //cout<<" Adding cell at "<<shiftx<<","<<shifty<<"\n";
                    idsNeighbour[0].push_back(id);
                    idsNeighbour[1].push_back(id);
                    idsNeighbour[2].push_back(id);
                    idsNeighbour[3].push_back(id);
                    continue;
                }

                if(r==0)
                {
                    DetId id = event().hgcalNavigator().topology()->offsetBy(detid, shiftx,  shifty);
                    const SimHit* hit   = (id!=DetId(0) ? event().simhit(id) : 0);

                    if     (shiftx==-n && shifty>-n ) { energy[0] += (hit ? hit->energy() : 0.); idsNeighbour[0].push_back(id); }
                    else if(shiftx<n   && shifty==n ) { energy[0] += (hit ? hit->energy() : 0.); idsNeighbour[0].push_back(id); }
                    //                                                                                                 
                    if     (shiftx>-n  && shifty==n ) { energy[1] += (hit ? hit->energy() : 0.); idsNeighbour[1].push_back(id); }
                    else if(shiftx==n  && shifty>-n ) { energy[1] += (hit ? hit->energy() : 0.); idsNeighbour[1].push_back(id); }
                    //                                                                                                 
                    if     (shiftx==n  && shifty< n ) { energy[2] += (hit ? hit->energy() : 0.); idsNeighbour[2].push_back(id); }
                    else if(shiftx>-n  && shifty==-n) { energy[2] += (hit ? hit->energy() : 0.); idsNeighbour[2].push_back(id); }
                    //                                                                                                 
                    if     (shiftx<n   && shifty==-n) { energy[3] += (hit ? hit->energy() : 0.); idsNeighbour[3].push_back(id); }
                    else if(shiftx==-n && shifty< n ) { energy[3] += (hit ? hit->energy() : 0.); idsNeighbour[3].push_back(id); }
                }
            }
        }
    }


    // find the max energy
    int maxIndex = 0;
    double maxEnergy = 0.;
    for(unsigned i=0;i<4;i++)
    {
        //cout<<" Energy "<<i<<" = "<<energy[i]<<"\n";
        if(energy[i]>maxEnergy)
        {
            maxEnergy = energy[i];
            maxIndex = i;
        }
    }


    // fill the 4 DetIds corresponding to the max 2x2
    vector<HGCEEDetId> ids = detids;
    for(auto itr=idsNeighbour[maxIndex].begin(); itr!=idsNeighbour[maxIndex].end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );


    return ids;

}

/*****************************************************************/
vector<HGCEEDetId> AnalysisSeeding::choose2x2( const vector<HGCEEDetId>& detids)
/*****************************************************************/
{
    vector<double> energy(4);
    energy[0] = 0.;
    energy[1] = 0.;
    energy[2] = 0.;
    energy[3] = 0.;
    vector<DetId> idsNW;
    vector<DetId> idsN;  
    vector<DetId> idsNE;
    vector<DetId> idsE ; 
    vector<DetId> idsSE ;
    vector<DetId> idsS  ;
    vector<DetId> idsSW ;
    vector<DetId> idsW  ;


    for(auto itr=detids.cbegin(); itr!=detids.end(); itr++)
    {
        const HGCEEDetId& detid = *itr;
        if(detid==DetId(0)) continue;

        // retrieve neighbour ids and hits
        DetId idNW = event().hgcalNavigator().topology()->offsetBy(detid, -1,  1);
        DetId idN  = event().hgcalNavigator().topology()->offsetBy(detid,  0,  1);
        DetId idNE = event().hgcalNavigator().topology()->offsetBy(detid,  1,  1);
        DetId idE  = event().hgcalNavigator().topology()->offsetBy(detid,  1,  0);
        DetId idSE = event().hgcalNavigator().topology()->offsetBy(detid,  1, -1);
        DetId idS  = event().hgcalNavigator().topology()->offsetBy(detid,  0, -1);
        DetId idSW = event().hgcalNavigator().topology()->offsetBy(detid, -1, -1);
        DetId idW  = event().hgcalNavigator().topology()->offsetBy(detid, -1,  0);

        idsNW.push_back(idNW);
        idsN .push_back(idN );
        idsNE.push_back(idNE);
        idsE .push_back(idE );
        idsSE.push_back(idSE);
        idsS .push_back(idS );
        idsSW.push_back(idSW);
        idsW .push_back(idW );

        //const SimHit* hitSeed = event().simhit(detid);
        const SimHit* hitNW   = (idNW!=DetId(0) ? event().simhit(idNW) : 0);
        const SimHit* hitN    = (idN !=DetId(0) ? event().simhit(idN ) : 0);
        const SimHit* hitNE   = (idNE!=DetId(0) ? event().simhit(idNE) : 0);
        const SimHit* hitE    = (idE !=DetId(0) ? event().simhit(idE ) : 0);
        const SimHit* hitSE   = (idSE!=DetId(0) ? event().simhit(idSE) : 0);
        const SimHit* hitS    = (idS !=DetId(0) ? event().simhit(idS ) : 0);
        const SimHit* hitSW   = (idSW!=DetId(0) ? event().simhit(idSW) : 0);
        const SimHit* hitW    = (idW !=DetId(0) ? event().simhit(idW ) : 0);

        // build the four 2x2 possible energies
        //energy[0] += (hitSeed ? hitSeed->energy() : 0.);
        energy[0] += (hitNW   ? hitNW->energy()   : 0.);
        energy[0] += (hitN    ? hitN ->energy()   : 0.);
        energy[0] += (hitW    ? hitW ->energy()   : 0.);
        //
        //energy[1] += (hitSeed ? hitSeed->energy() : 0.);
        energy[1] += (hitNE   ? hitNE->energy()   : 0.);
        energy[1] += (hitN    ? hitN ->energy()   : 0.);
        energy[1] += (hitE    ? hitE ->energy()   : 0.);
        //
        //energy[2] += (hitSeed ? hitSeed->energy() : 0.);
        energy[2] += (hitSE   ? hitSE->energy()   : 0.);
        energy[2] += (hitS    ? hitS ->energy()   : 0.);
        energy[2] += (hitE    ? hitE ->energy()   : 0.);
            //
        //energy[3] += (hitSeed ? hitSeed->energy() : 0.);
        energy[3] += (hitSW   ? hitSW->energy()   : 0.);
        energy[3] += (hitS    ? hitS ->energy()   : 0.);
        energy[3] += (hitW    ? hitW ->energy()   : 0.);


    }


    // find the max energy
    int maxIndex = -1;
    double maxEnergy = 0.;
    for(unsigned i=0;i<4;i++)
    {
        cout<<" Energy "<<i<<" = "<<energy[i]<<"\n";
        if(energy[i]>maxEnergy)
        {
            maxEnergy = energy[i];
            maxIndex = i;
        }
    }


    // fill the 4 DetIds corresponding to the max 2x2
    vector<HGCEEDetId> ids = detids;
    if(maxIndex==0)
    {
        for(auto itr=idsNW.begin(); itr!=idsNW.end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsN .begin(); itr!=idsN .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsW .begin(); itr!=idsW .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
    }
    else if(maxIndex==1)
    {
        for(auto itr=idsNE.begin(); itr!=idsNE.end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsN .begin(); itr!=idsN .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsE .begin(); itr!=idsE .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
    }
    else if(maxIndex==2)
    {
        for(auto itr=idsSE.begin(); itr!=idsSE.end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsS .begin(); itr!=idsS .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsE .begin(); itr!=idsE .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
    }
    else if(maxIndex==3)
    {
        for(auto itr=idsSW.begin(); itr!=idsSW.end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsS .begin(); itr!=idsS .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
        for(auto itr=idsW .begin(); itr!=idsW .end(); itr++) ids.push_back( (HGCEEDetId)(*itr) );
    }




    return ids;

}



/*****************************************************************/
void AnalysisSeeding::findMaxHit(int i)
/*****************************************************************/
{
    m_chosenHit = NULL;

    if(m_sample=="electron")
    {
        m_genParticle = &event().genparticles()[i];
        double eta = m_genParticle->Eta();
        double phi = m_genParticle->Phi();
        double maxE = 0.;
        for(auto itrHit=event().simhits().cbegin(); itrHit!=event().simhits().end(); itrHit++)
        {
            const SimHit& hit = *itrHit;
            // match only to hits in layers 10-20
            if(hit.layer()<10 || hit.layer()>20) continue;

            double deta = (hit.eta() - eta);
            double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(dr<0.3 && hit.energy()>maxE)
            {
                m_chosenHit = &hit;
                maxE = hit.energy();
            }
        }
    }
    else if(m_sample=="electronPU")
    {
        m_genParticle = &event().genparticles()[i];
        double eta = m_genParticle->Eta();
        double phi = m_genParticle->Phi();
        double maxE = 0.;
        for(auto itrHit=event().hardsimhits().cbegin(); itrHit!=event().hardsimhits().end(); itrHit++)
        {
            const SimHit& hit = *itrHit;
            // match only to hits in layers 10-20
            if(hit.layer()<10 || hit.layer()>20) continue;

            double deta = (hit.eta() - eta);
            double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
            double dr = sqrt(deta*deta + dphi*dphi);
            if(dr<0.1 && hit.energy()>maxE)
            {
                m_chosenHit = &hit;
                maxE = hit.energy();
            }
        }
    }
    else // MinBias. Choose random direction
    {
        //while(!m_chosenHit)
        //{
        //    double eta = m_random.Uniform(1.5, 3.);
        //    double phi = m_random.Uniform(-M_PI, M_PI);
        //    double maxE = 0.;
        //    for(auto itrHit=event().simhits().cbegin(); itrHit!=event().simhits().end(); itrHit++)
        //    {
        //        const SimHit& hit = *itrHit;
        //        double deta = (hit.eta() - eta);
        //        double dphi = TVector2::Phi_mpi_pi(hit.phi() - phi);
        //        double dr = sqrt(deta*deta + dphi*dphi);
        //        if(dr<0.2 && hit.energy()>maxE)
        //        {
        //            m_chosenHit = &hit;
        //            maxE = hit.energy();
        //        }
        //    }
        //}
        if((unsigned)i>=event().sortedHits().size()) return;
        m_chosenHit = event().sortedHits()[i];
        // make sure a close-by hit has not been already picked
        for(auto itr=m_usedHits.cbegin(); itr!=m_usedHits.end();itr++)
        {
            const SimHit* usedHit = *itr;
            double deta = (m_chosenHit->eta() - usedHit->eta());
            double dphi = TVector2::Phi_mpi_pi(m_chosenHit->phi() - usedHit->phi());
            double dr = sqrt(deta*deta + dphi*dphi);
            if(dr<0.1) 
            {
                m_chosenHit = 0;
                return;
            }
        }
        m_usedHits.push_back(m_chosenHit);

        //cout<<" E="<<m_chosenHit->energy()<<", eta="<<m_chosenHit->eta()<<", phi="<<m_chosenHit->phi()<<",layer="<<m_chosenHit->layer()<<",cell="<<m_chosenHit->cell()<<"\n";
    }

}

