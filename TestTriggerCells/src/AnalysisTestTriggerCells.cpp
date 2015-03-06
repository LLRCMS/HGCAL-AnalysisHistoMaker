/**
 *  @file  AnalysisTestTriggerCells.cpp
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    02/03/2015
 *
 *  @internal
 *     Created :  02/03/2015
 * Last update :  02/03/2015 11:22:06
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */




#include <iostream>
#include <fstream>
#include <stdexcept>

#include "TVector2.h"

#include "AnHiMaHGCAL/TestTriggerCells/interface/AnalysisTestTriggerCells.h"
#include "AnHiMaHGCAL/Core/interface/Utilities.h"
#include "AnHiMaHGCAL/HGCALCommon/interface/Constants.h"





using namespace AnHiMa;
using namespace std;



/*****************************************************************/
AnalysisTestTriggerCells::AnalysisTestTriggerCells():IAnalysis()
/*****************************************************************/
{

}



/*****************************************************************/
AnalysisTestTriggerCells::~AnalysisTestTriggerCells()
/*****************************************************************/
{

}



/*****************************************************************/
bool AnalysisTestTriggerCells::initialize(const string& parameterFile)
/*****************************************************************/
{
    bool status = IAnalysis::initialize(parameterFile);
    if(!status) return status;
    event().connectVariables(m_inputChain);

    m_cellToTriggerCell.resize(30);


    // Read trigger cell mapping
    ifstream ifs("/home/llr/cms/sauvan/CMSSW/HGCAL/CMSSW_6_2_0_SLHC20/src/AnHiMaHGCAL/HGCALCommon/data/cellsToTriggerCellsMap.txt", std::ifstream::in);
    if(!ifs.is_open()) 
    {
        cerr<<"ERROR: Cannot open trigger cell mapping \n";
        return false;
    }
    short layer = 0;
    short cell = 0;
    short triggerCell = 0;
    short subsector = 0;
    
    for(; ifs>>layer>>cell>>triggerCell>>subsector; )
    {
        if(layer>=30 || layer<0) 
        {
            cerr<<"ERROR: TriggerHitFactory::initialize(): Found layer '"<<layer<<"' in trigger cell mapping. Discard it.\n";
            continue;
        }
        short newcell = (subsector==1 ? cell+1 : -(cell+1));
        auto ret = m_cellToTriggerCell[layer].insert( std::make_pair(newcell, triggerCell) );
        if(ret.second==false)
        {
            cerr<<"ERROR: TriggerHitFactory::initialize(): Duplicated entry ("<<layer<<","<<newcell<<") in trigger cell mapping. Use first mapping found.\n";
        }
        vector<short> vec(0);
        m_triggerCellToCell[layer].insert( std::make_pair(triggerCell, vec) );
        m_triggerCellToCell[layer].at(triggerCell).push_back(newcell);
    } 
    assert(ifs.eof());
    ifs.close();

    // compute positions of trigger cells
    for(unsigned sec=1;sec<=18; sec++)
    {
        for(unsigned layer=1; layer<=30; layer++)
        {
            for(const auto& trigger_cells : m_triggerCellToCell[layer-1])
            {
                const short triggerCell = trigger_cells.first;
                const vector<short>& cells = trigger_cells.second;
                SimplePosition pos;
                pos.x   = 0.;
                pos.y   = 0.;
                for(const auto& c : cells)
                {
                    int subsec = (c>0 ? 1 : -1);
                    int cell = abs(c)-1;
                    HGCEEDetId detid(HGCEE, 1, layer, sec, subsec, cell);
                    const GlobalPoint point = event().hgcalNavigator().geometry()->getPosition(detid);
                    //pos.eta += point.eta();
                    //pos.phi += (pos.phi + TVector2::Phi_mpi_pi(point.phi() - pos.phi));
                    pos.x   += point.x();
                    pos.y   += point.y();
                }
                //pos.eta /= (float)cells.size();
                //pos.phi  = TVector2::Phi_mpi_pi(pos.phi/(float)cells.size());
                pos.x   /= (float)cells.size();
                pos.y   /= (float)cells.size();
                m_triggerCellPositions[sec-1][layer-1].insert( std::make_pair(triggerCell, pos) );
            }
        }
    }


    return true;
}


/*****************************************************************/
void AnalysisTestTriggerCells::execute()
/*****************************************************************/
{
    cerr<<"Event "<<event().event()<<"\n";
    event().update();
    if(!event().passSelection()) return;


    testTriggerCellBuilding();
    testTriggerCellNavigation();

}


/*****************************************************************/
void AnalysisTestTriggerCells::testTriggerCellBuilding()
/*****************************************************************/
{
    // create trigger cells as Tower objects
    map<int, Tower> triggerCells;
    for(const auto& hit : event().simhits() )
    {
        short layer = hit.layer()-1;
        short cell = (hit.cell()+1)*hit.subsector();
        short triggerCell = m_cellToTriggerCell[layer].at(cell);
        int index = cellIndex(hit.zside(), hit.layer(), hit.sector(), triggerCell);
        Tower tower;
        triggerCells.insert( pair<int, Tower>( index, tower ));
        triggerCells.at(index).addHit(hit);
    }

    if(event().triggerhits().size()!=triggerCells.size())
    {
        cout<<"Problem: inconsistent number of trigger hits "<<event().triggerhits().size()<<"!="<<triggerCells.size()<<"\n";
    }
    // compare trigger hits with previous towers
    for(const auto& hit : event().triggerhits())
    {
        int index = cellIndex(hit.zside(), hit.layer(), hit.sector(), hit.cell());
        try
        {
            const Tower& tower = triggerCells.at(index);
            if(hit.energy()!=tower.energy())
            {
                cout<<"Problem in trigger cell layer="<<hit.layer()<<", sector="<<hit.sector()<<", cell="<<hit.cell()<<", energy="<<hit.energy()<<"!="<<tower.energy()<<"\n";
            }
            double towerEta = tower.eta();
            double towerPhi = tower.phi();
            double towerX   = tower.x();
            double towerY   = tower.y();
            double triggerEta = hit.eta();
            double triggerPhi = hit.phi();
            double triggerX   = hit.x();
            double triggerY   = hit.y();
            unsigned nCells = m_triggerCellToCell[hit.layer()-1].at(hit.cell()).size();
            if( fabs(triggerEta-towerEta)>0.02 ) 
            {
                cout<<"Eta difference: eta: "<<triggerEta<<" "<<towerEta<<" "<<triggerEta-towerEta;
                cout<<" ncells="<<nCells<<"\n";
                cout<<"  trigger cell "<<hit.cell()<<" sector "<<hit.sector()<<" layer "<<hit.layer()<<"\n";
                cout<<"  cells eta = {";
                for(const auto& c : m_triggerCellToCell[hit.layer()-1].at(hit.cell()))
                {
                    int subsec = (c>0 ? 1 : -1);
                    int cell = abs(c)-1;
                    HGCEEDetId detid(HGCEE, hit.zside(), hit.layer(), hit.sector(), subsec, cell);
                    const GlobalPoint point = event().hgcalNavigator().geometry()->getPosition(detid); 
                    cout<<point.eta()<<" ";
                }
                cout<<"}\n";
                
            }
            if( (fabs(TVector2::Phi_mpi_pi(triggerPhi-towerPhi))>0.02 && nCells<=4) || (fabs(TVector2::Phi_mpi_pi(triggerPhi-towerPhi))>0.03 && nCells>4))
            {
                cout<<"Phi difference: phi: "<<triggerPhi<<" "<<towerPhi<<" "<<triggerPhi-towerPhi;
                cout<<" ncells="<<nCells<<"\n";
                cout<<"  trigger cell "<<hit.cell()<<" sector "<<hit.sector()<<" layer "<<hit.layer()<<"\n";
                cout<<"  cells phi = {";
                for(const auto& c : m_triggerCellToCell[hit.layer()-1].at(hit.cell()))
                {
                    int subsec = (c>0 ? 1 : -1);
                    int cell = abs(c)-1;
                    HGCEEDetId detid(HGCEE, hit.zside(), hit.layer(), hit.sector(), subsec, cell);
                    const GlobalPoint point = event().hgcalNavigator().geometry()->getPosition(detid); 
                    cout<<point.phi()<<" ";
                }
                cout<<"}\n";
            }
            if( (fabs(triggerX-towerX)>1. && nCells<=4) || (fabs(triggerX-towerX)>2. && nCells>5) )
            {
                cout<<"X difference:  x: "<<triggerX<<" "<<towerX<<" "<<triggerX-towerX;
                cout<<" ncells="<<nCells<<"\n";
                cout<<"  trigger cell "<<hit.cell()<<" sector "<<hit.sector()<<" layer "<<hit.layer()<<"\n";
                cout<<"  cells X = {";
                for(const auto& c : m_triggerCellToCell[hit.layer()-1].at(hit.cell()))
                {
                    int subsec = (c>0 ? 1 : -1);
                    int cell = abs(c)-1;
                    HGCEEDetId detid(HGCEE, hit.zside(), hit.layer(), hit.sector(), subsec, cell);
                    const GlobalPoint point = event().hgcalNavigator().geometry()->getPosition(detid); 
                    cout<<point.x()<<" ";
                }
                cout<<"}\n";
            }
            if( (fabs(triggerY-towerY)>1 && nCells<=4) || (fabs(triggerY-towerY)>2. && nCells>5))
            {
                cout<<"Y difference:  y: "<<triggerY<<" "<<towerY<<" "<<triggerY-towerY;
                cout<<" ncells="<<nCells<<"\n";
                cout<<"  trigger cell "<<hit.cell()<<" sector "<<hit.sector()<<" layer "<<hit.layer()<<"\n";
                cout<<"  cells Y = {";
                for(const auto& c : m_triggerCellToCell[hit.layer()-1].at(hit.cell()))
                {
                    int subsec = (c>0 ? 1 : -1);
                    int cell = abs(c)-1;
                    HGCEEDetId detid(HGCEE, hit.zside(), hit.layer(), hit.sector(), subsec, cell);
                    const GlobalPoint point = event().hgcalNavigator().geometry()->getPosition(detid); 
                    cout<<point.y()<<" ";
                }
                cout<<"}\n";
            }

        } catch(std::out_of_range& e)
        {
            cout<<"Problem: inconsistent set of trigger cells: cannot find cell layer="<<hit.layer()<<", sector="<<hit.sector()<<", cell="<<hit.cell()<<"\n";
        }

    }

    // Look at cell distances inside trigger cells
    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    for(const auto& hit : event().triggerhits())
    {
        double maxdxdy = 0.;
        double dxmax = 0.;
        double dymax = 0.;
        const auto& cellList =  m_triggerCellToCell[hit.layer()-1].at(hit.cell());
        for(unsigned i=0;i<cellList.size()-1; i++)
        {
            int c1 = cellList[i];
            int subsec1 = (c1>0 ? 1 : -1);
            int cell1 = abs(c1)-1;
            HGCEEDetId detid1(HGCEE, hit.zside(), hit.layer(), hit.sector(), subsec1, cell1);
            const GlobalPoint point1 = event().hgcalNavigator().geometry()->getPosition(detid1); 
            double x1 = point1.x();
            double y1 = point1.y();
            for(unsigned j=i;j<cellList.size(); j++)
            {
                int c2 = cellList[j];
                int subsec2 = (c2>0 ? 1 : -1);
                int cell2 = abs(c2)-1;
                HGCEEDetId detid2(HGCEE, hit.zside(), hit.layer(), hit.sector(), subsec2, cell2);
                const GlobalPoint point2 = event().hgcalNavigator().geometry()->getPosition(detid2); 
                double x2 = point2.x();
                double y2 = point2.y();
                //
                double dx = (x2-x1);
                double dy = (y2-y1);
                if(sqrt(dx*dx+dy*dy)>maxdxdy)
                {
                    maxdxdy = sqrt(dx*dx+dy*dy);
                    dxmax = dx;
                    dymax = dy;
                }
            }
            if(maxdxdy>10.)
            {
                cout<<"trigger cell "<<hit.cell()<<" layer "<<hit.layer()<<" max DR(cells)="<<maxdxdy<<"\n";
            }
            m_histos.Fill1BinHisto(10+hoffset, cellList.size(), fabs(dxmax), fabs(dymax), weight, sysNum);
        }
    }

}

/*****************************************************************/
void AnalysisTestTriggerCells::testTriggerCellNavigation()
/*****************************************************************/
{
    short sysNum = 0;
    float weight = 1.;
    int hoffset  = 0;

    for(const auto& hit : event().triggerhits())
    {
        if(hit.sector()!=1) continue; // Only sector parallel to the x-axis
        HGCEEDetId id(HGCEE, hit.zside(), hit.layer(), hit.sector(), hit.subsector(), hit.cell());
        const auto neighborsEast  = event().hgcalNavigator().eastTrigger(id);
        const auto neighborsWest  = event().hgcalNavigator().westTrigger(id);
        const auto neighborsNorth = event().hgcalNavigator().northTrigger(id);
        const auto neighborsSouth = event().hgcalNavigator().southTrigger(id);

        double x1 = hit.x();
        double y1 = hit.y();
        const auto& cellList =  m_triggerCellToCell[hit.layer()-1].at(hit.cell());
        //
        for(const auto& neigh : neighborsEast)
        {
            cerr<<neigh.sector()<<" "<<neigh.layer()<<" "<<neigh.cell()<<"\n";
            SimplePosition pos = m_triggerCellPositions.at(neigh.sector()-1).at(neigh.layer()-1).at(neigh.cell());
            double x2 = pos.x;
            double y2 = pos.y;
            //
            double dx = (x2-x1);
            double dy = (y2-y1);
            m_histos.Fill1BinHisto(20+hoffset, cellList.size(), dx, dy, weight, sysNum);
        }
        //
        for(const auto& neigh : neighborsWest)
        {
            cerr<<neigh.sector()<<" "<<neigh.layer()<<" "<<neigh.cell()<<"\n";
            SimplePosition pos = m_triggerCellPositions.at(neigh.sector()-1).at(neigh.layer()-1).at(neigh.cell());
            double x2 = pos.x;
            double y2 = pos.y;
            //
            double dx = (x2-x1);
            double dy = (y2-y1);
            m_histos.Fill1BinHisto(30+hoffset, cellList.size(), dx, dy, weight, sysNum);
        }
        //
        for(const auto& neigh : neighborsNorth)
        {
            cerr<<neigh.sector()<<" "<<neigh.layer()<<" "<<neigh.cell()<<"\n";
            SimplePosition pos = m_triggerCellPositions.at(neigh.sector()-1).at(neigh.layer()-1).at(neigh.cell());
            double x2 = pos.x;
            double y2 = pos.y;
            //
            double dx = (x2-x1);
            double dy = (y2-y1);
            m_histos.Fill1BinHisto(40+hoffset, cellList.size(), dx, dy, weight, sysNum);
        }
        //
        for(const auto& neigh : neighborsSouth)
        {
            cerr<<neigh.sector()<<" "<<neigh.layer()<<" "<<neigh.cell()<<"\n";
            SimplePosition pos = m_triggerCellPositions.at(neigh.sector()-1).at(neigh.layer()-1).at(neigh.cell());
            double x2 = pos.x;
            double y2 = pos.y;
            //
            double dx = (x2-x1);
            double dy = (y2-y1);
            m_histos.Fill1BinHisto(50+hoffset, cellList.size(), dx, dy, weight, sysNum);
        }







        //cout<<"sector "<<hit.sector()<<" layer "<<hit.layer()<<" cell "<<hit.cell()<<"\n";
        //cout<<" Neighbors east = ";
        //for(const auto& neigh : neighborsEast)
        //{
        //    cout<<"("<<neigh.sector()<<","<<neigh.layer()<<","<<neigh.cell()<<")";
        //}
        //cout<<"\n";
        //cout<<" Neighbors west = ";
        //for(const auto& neigh : neighborsWest)
        //{
        //    cout<<"("<<neigh.sector()<<","<<neigh.layer()<<","<<neigh.cell()<<")";
        //}
        //cout<<"\n";
        //cout<<" Neighbors north = ";
        //for(const auto& neigh : neighborsNorth)
        //{
        //    cout<<"("<<neigh.sector()<<","<<neigh.layer()<<","<<neigh.cell()<<")";
        //}
        //cout<<"\n";
        //cout<<" Neighbors south = ";
        //for(const auto& neigh : neighborsSouth)
        //{
        //    cout<<"("<<neigh.sector()<<","<<neigh.layer()<<","<<neigh.cell()<<")";
        //}
        //cout<<"\n";
    }

}



/*****************************************************************/
unsigned AnalysisTestTriggerCells::cellIndex(int zside, int layer, int sector, int cell) const
/*****************************************************************/
{
    int z = (zside==-1 ? 0 : 1);
    int l = layer-1;
    int s = sector-1;
    //int ss = (subsector==-1 ? 0 : 1);
    int c = cell;
    //return (unsigned)(c + ss*3600 + s*7200 + l*129600 + z*3888000);
    return (unsigned)(c + s*3600 + l*64800 + z*1944000);
}


