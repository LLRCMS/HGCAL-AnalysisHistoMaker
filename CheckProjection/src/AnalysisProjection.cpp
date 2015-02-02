/**
 *  @file  AnalysisProjection.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    23/01/2015
 *
 *  @internal
 *     Created :  23/01/2015
 * Last update :  23/01/2015 12:24:47
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */





#include <iostream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/CheckProjection/interface/AnalysisProjection.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"







using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisProjection::AnalysisProjection():IAnalysis(),
    m_eventCount(0)
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisProjection::~AnalysisProjection()
/*****************************************************************/
{

}



/*****************************************************************/
bool AnalysisProjection::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    // will built projected hits onto the HGCEE layer 1
    string projectionMappingFile(m_reader.params().GetValue("ProjectionMapping", ""));
    event().initializeHitProjection(projectionMappingFile);
    fillProjectionMapping(projectionMappingFile);


    return true;
}


/*****************************************************************/
void AnalysisProjection::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;


    fillHistos();

    m_eventCount++;

}

/*****************************************************************/
void AnalysisProjection::fillHistos()
/*****************************************************************/
{

    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    for(const auto& hit : event().simhits())
    {
        Tower tower;
        TowerCalibrator calib;
        tower.addHit(hit);
        calib.calibrate(tower);
        m_histos.Fill1BinHisto(100+hoffset, m_eventCount, hit.x(), hit.y(), weight*tower.calibratedEnergy(), sysNum);
    }

    for(const auto& projhit : event().projectedhits())
    {
        m_histos.Fill1BinHisto(200+hoffset, m_eventCount, projhit.x(), projhit.y(), weight*projhit.energy(), sysNum);
        // compare energies around high energy hits
        if(projhit.energy()>5.)
        {
            // find position of projected hit
            double xproj = projhit.x();
            double yproj = projhit.y();
            double zsideproj = projhit.zside();
            // build tower around the projected hit
            Tower tower;
            TowerCalibrator calib;
            for(const auto& hit : event().simhits())
            {
                double dx = fabs(xproj - hit.x());
                double dy = fabs(yproj - hit.y());

                if(dx<0.475 && dy<0.475 && hit.zside()==zsideproj)
                {
                    tower.addHit(hit);
                }
            }
            calib.calibrate(tower);

            double reducedPhi = (projhit.phi()+TMath::Pi()/18.)/(2*TMath::Pi()/18.); // divide by 20°
            double localPhi = (reducedPhi - floor(reducedPhi))*(2*TMath::Pi()/18.);
            m_histos.Fill1BinHisto(10+hoffset, fabs(projhit.eta()), projhit.energy()/tower.calibratedEnergy(), weight, sysNum);
            m_histos.Fill1BinHisto(20+hoffset, localPhi, projhit.energy()/tower.calibratedEnergy(), weight, sysNum);
            
            /// Towers with more energy than projected hits
            if(projhit.energy()/tower.calibratedEnergy() < 1.)
            {
                m_histos.FillHisto(300+hoffset, projhit.x(), projhit.y(), weight, sysNum);
            }
            /// Towers with less energy than projected hits
            if(projhit.energy()/tower.calibratedEnergy() > 1.)
            {
                m_histos.FillHisto(301+hoffset, projhit.x(), projhit.y(), weight, sysNum);
            }

        }
    }

}

/*****************************************************************/
HGCEEDetId AnalysisProjection::projectCell(const HGCEEDetId& cellid)
/*****************************************************************/
{
    int layer = cellid.layer();
    int cell  = cellid.cell();
    std::pair<unsigned,unsigned> layer_cell = make_pair(layer,cell);
    const auto& itr = m_projectionMapping.find(layer_cell);
    if(itr==m_projectionMapping.end())
    {
        cerr<<"WARNING: Cannot find projected cell for (layer,cell)=("<<layer<<","<<cell<<")\n";
        return DetId(0);
    }
    return HGCEEDetId(HGCEE, cellid.zside(), 1, cellid.sector(), cellid.subsector(), itr->second);
}


/*****************************************************************/
void AnalysisProjection::fillProjectionMapping(const string& inputFile)
/*****************************************************************/
{
    fstream input(inputFile,fstream::in);
    if(!input.is_open())
    {
        cerr<<"ERROR: Cannot open projection mapping file '"<<inputFile<<"'\n";
        return;
    }
    while(!input.eof())
    {
        int layer, cell, projcell;
        input >> layer >> cell >> projcell;
        m_projectionMapping[make_pair(layer,cell)] = projcell;
    }
    input.close();

}






