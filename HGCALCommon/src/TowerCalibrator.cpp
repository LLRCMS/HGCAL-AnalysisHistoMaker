


#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

#include <math.h>
#include <iostream>



using namespace std;
using namespace AnHiMa;


/*****************************************************************/
TowerCalibrator::TowerCalibrator():
    // EE
    m_mipValueInGeV_EE(55.1*1e-6),
    m_coeff_a_EE(82.8),
    m_coeff_b_EE(1e6),
    m_coeff_c_EE(1e6),
    m_weights_EE(30),
    // HEF
    m_mipValueInGeV_HEF(85.0*1e-6),
    m_coeff_a_HEF(1.0),
    m_coeff_b_HEF(1e6),
    m_coeff_c_HEF(1e6),
    m_weights_HEF(12),
    // HEB
    m_mipValueInGeV_HEB(1498.4*1e-6),
    m_coeff_a_HEB(1.0),
    m_coeff_b_HEB(1e6),
    m_coeff_c_HEB(1e6),
    m_weights_HEB(12),
    // 1/2 layer calibration
    m_halfLayersCalib(15)
/*****************************************************************/
{
    // EE
    for(unsigned l=0; l<m_weights_EE.size(); l++)
    {
        double weight = 1.;
        if(l==0) weight = 0.080;
        else if(l==1) weight = 0.92;
        else if(l>=2 && l<=10) weight = 0.62;
        else if(l>=11 && l<=20) weight = 0.81;
        else if(l>=21 && l<=29) weight = 1.19;
        m_weights_EE[l] = weight;
    }
    // HEF
    for(unsigned l=0; l<m_weights_HEF.size(); l++)
    {
        double weight = 1.;
        if(l==0) weight = 0.0464;
        else if(l>=1 && l<=11) weight = 0.0474;
        m_weights_HEF[l] = weight;
    }
    // HEB
    for(unsigned l=0; l<m_weights_HEB.size(); l++)
    {
        double weight = 1.;
        if(l<=11) weight = 0.1215;
        m_weights_HEB[l] = weight;
    }
    // 1/2 layer calibration
    // /home/sauvan/Documents/postdoc_LLR/ANALYSES/HGCAL/Plots/LayerCalibration/correctionsHisto/layerCalibrations.txt 
    m_halfLayersCalib[0]  = 1.02071661478; // 2 
    m_halfLayersCalib[1]  = 1.62790603275; // 4
    m_halfLayersCalib[2]  = 1.74677030342; // 6
    m_halfLayersCalib[3]  = 1.83534112526; // 8
    m_halfLayersCalib[4]  = 1.90645895385; // 10
    m_halfLayersCalib[5]  = 1.71465184791; // 12
    m_halfLayersCalib[6]  = 2.03674910703; // 14
    m_halfLayersCalib[7]  = 2.10002662805; // 16
    m_halfLayersCalib[8]  = 2.14546184065; // 18
    m_halfLayersCalib[9]  = 2.17769377047; // 20
    m_halfLayersCalib[10] = 1.88632355541; // 22
    m_halfLayersCalib[11] = 2.37351161737; // 24
    m_halfLayersCalib[12] = 2.3913787003;  // 26
    m_halfLayersCalib[13] = 2.4028844731;  // 28
    m_halfLayersCalib[14] = 2.42617555279; // 30
}



/*****************************************************************/
TowerCalibrator::~TowerCalibrator()
/*****************************************************************/
{
}



/*****************************************************************/
void TowerCalibrator::calibrate(Tower& tower, bool halflayers)
/*****************************************************************/
{
    const vector<const SimHit*>& hits = tower.hits();
    // Taken from  https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/RecoParticleFlow/PFClusterProducer/src/HGCEEElectronEnergyCalibrator.cc
    // and https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/RecoParticleFlow/PFClusterProducer/src/HGCHEHadronicEnergyCalibrator.cc
    double eCorr = 0.0;
    vector< vector<double> > eCorrLayer;
    eCorrLayer.push_back(vector<double>(31)); // EE
    eCorrLayer.push_back(vector<double>(13)); // HEF
    eCorrLayer.push_back(vector<double>(13)); // HEB
    for(unsigned s=0;s<eCorrLayer.size();s++)
    {
        for(unsigned l=0;l<eCorrLayer[s].size();l++)
        {
            eCorrLayer[s][l] = 0.;
        }
    }
    double eta = tower.eta();
    for( const auto& h : hits ) 
    {
        const SimHit& hit = *h;
        int layer = hit.layer();
        int subdet = hit.subdet();
        const double effMIP_to_InvGeV_EE  = m_coeff_a_EE /(1.0 + exp(-m_coeff_c_EE  - m_coeff_b_EE *cosh(eta)));
        const double effMIP_to_InvGeV_HEF = m_coeff_a_HEF/(1.0 + exp(-m_coeff_c_HEF - m_coeff_b_HEF*cosh(eta)));
        const double effMIP_to_InvGeV_HEB = m_coeff_a_HEB/(1.0 + exp(-m_coeff_c_HEB - m_coeff_b_HEB*cosh(eta))); 
        double mip_value = 0.0;
        double mip2gev = 0.0;
        double weight = 0.0;
        switch(subdet)
        {
            case 3: // EE
                mip_value = m_mipValueInGeV_EE;
                mip2gev = effMIP_to_InvGeV_EE;
                weight = m_weights_EE[layer-1];
                break;
            case 4: // HEF
                mip_value = m_mipValueInGeV_HEF;
                mip2gev = effMIP_to_InvGeV_HEF;
                weight = m_weights_HEF[layer-1];
                break;
            case 5: // HEB
                mip_value = m_mipValueInGeV_HEB;
                mip2gev = effMIP_to_InvGeV_HEB;
                weight = m_weights_HEB[layer-1];
                break;
            default:
                break;
        }
        const double energy_MIP = hit.energy()/mip_value; 
        if(halflayers && layer%2!=0) cout<<"WARNING: 1/2 layer calibration requested but found an odd layer ("<<layer<<")\n";
        const double corr = (halflayers ? m_halfLayersCalib[(layer-1)/2] : 1.);
        eCorr += weight*energy_MIP/mip2gev*corr;
        eCorrLayer[subdet-3][layer] += weight*energy_MIP/mip2gev*corr;

    }



    tower.setCalibratedEnergy(eCorr);
    for(unsigned s=0;s<eCorrLayer.size();s++)
    {
        for(unsigned l=0;l<eCorrLayer[s].size();l++)
        {
            tower.setLayerCalibratedEnergy(l, eCorrLayer[s][l], s+3);
        }
    }
}


/*****************************************************************/
double TowerCalibrator::calibratedEnergy(double energy, double eta, int layer, int subdet)
/*****************************************************************/
{
    const double effMIP_to_InvGeV_EE  = m_coeff_a_EE /(1.0 + exp(-m_coeff_c_EE  - m_coeff_b_EE *cosh(eta)));
    const double effMIP_to_InvGeV_HEF = m_coeff_a_HEF/(1.0 + exp(-m_coeff_c_HEF - m_coeff_b_HEF*cosh(eta)));
    const double effMIP_to_InvGeV_HEB = m_coeff_a_HEB/(1.0 + exp(-m_coeff_c_HEB - m_coeff_b_HEB*cosh(eta))); 
    double mip_value = 0.0;
    double mip2gev = 0.0;
    double weight = 0.0;
    switch(subdet)
    {
        case 3: // EE
            mip_value = m_mipValueInGeV_EE;
            mip2gev = effMIP_to_InvGeV_EE;
            weight = m_weights_EE[layer-1];
            break;
        case 4: // HEF
            mip_value = m_mipValueInGeV_HEF;
            mip2gev = effMIP_to_InvGeV_HEF;
            weight = m_weights_HEF[layer-1];
            break;
        case 5: // HEB
            mip_value = m_mipValueInGeV_HEB;
            mip2gev = effMIP_to_InvGeV_HEB;
            weight = m_weights_HEB[layer-1];
            break;
        default:
            break;
    }
    const double energy_MIP = energy/mip_value; 
    double eCorr = weight*energy_MIP/mip2gev;

    return eCorr;
}
