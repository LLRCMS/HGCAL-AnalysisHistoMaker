


#include "AnHiMaHGCAL/HGCALCommon/interface/TowerCalibrator.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/SimHit.h"

#include <math.h>



using namespace std;
using namespace AnHiMa;


/*****************************************************************/
TowerCalibrator::TowerCalibrator():
    m_mipValueInGeV(55.1*1e-6),
    m_coeff_a(82.8),
    m_coeff_b(1e6),
    m_coeff_c(1e6),
    m_weights(30)
/*****************************************************************/
{
    for(unsigned l=0; l<m_weights.size(); l++)
    {
        double weight = 1.;
        if(l==0) weight = 0.080;
        else if(l==1) weight = 0.92;
        else if(l>=2 && l<=10) weight = 0.62;
        else if(l>=11 && l<=20) weight = 0.81;
        else if(l>=21 && l<=29) weight = 1.19;
        m_weights[l] = weight;
    }
}



/*****************************************************************/
TowerCalibrator::~TowerCalibrator()
/*****************************************************************/
{
}



/*****************************************************************/
void TowerCalibrator::calibrate(Tower& tower)
/*****************************************************************/
{
    const vector<const SimHit*>& hits = tower.hits();
    // Taken from  https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/RecoParticleFlow/PFClusterProducer/src/HGCEEElectronEnergyCalibrator.cc
    double eCorr = 0.0;
    vector<double> eCorrLayer(31);
    for(unsigned l=0;l<eCorrLayer.size();l++)
    {
        eCorrLayer[l] = 0.;
    }
    double eta = tower.eta();
    for( const auto& h : hits ) 
    {
        const SimHit& hit = *h;
        int layer = hit.layer();
        double energy_MIP = hit.energy()/m_mipValueInGeV;
        eCorr += m_weights[layer-1]*energy_MIP;
        eCorrLayer[layer] += m_weights[layer-1]*energy_MIP;
    }

    double effMIP_to_InvGeV = m_coeff_a/(1.0 + exp(-m_coeff_c - m_coeff_b*cosh(eta)));

    tower.setCalibratedEnergy(eCorr/effMIP_to_InvGeV);
    for(unsigned l=1;l<eCorrLayer.size();l++)
    {
        tower.setLayerCalibratedEnergy(l, eCorrLayer[l]/effMIP_to_InvGeV);
    }
}


/*****************************************************************/
double TowerCalibrator::calibratedEnergy(double energy, double eta, int layer)
/*****************************************************************/
{
;
    double energy_MIP = energy/m_mipValueInGeV;
    double eCorr = m_weights[layer-1]*energy_MIP;
    double effMIP_to_InvGeV = m_coeff_a/(1.0 + exp(-m_coeff_c - m_coeff_b*cosh(eta)));

    return eCorr/effMIP_to_InvGeV;
}
