

#include <sstream>
#include <fstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/ProjectionMapping/interface/ProjectionMapProducer.h"

using namespace std;
using namespace AnHiMa;


const int ProjectionMapProducer::MAXCELLS = 2386;

/*****************************************************************/
ProjectionMapProducer::ProjectionMapProducer():
    m_output(0),
    m_projectedBins(0),
    m_unProjected(0),
    m_binToCell(MAXCELLS),
    m_xN(0),
    m_yN(0),
    m_xMin(0.),
    m_xMax(0.),
    m_yMin(0.),
    m_yMax(0.),
    m_cellsInLayers(30)
/*****************************************************************/
{
    m_cellsInLayers[0]  = 2117;
    m_cellsInLayers[1]  = 2124;
    m_cellsInLayers[2]  = 2124;
    m_cellsInLayers[3]  = 2126;
    m_cellsInLayers[4]  = 2126;
    m_cellsInLayers[5]  = 2159;
    m_cellsInLayers[6]  = 2159;
    m_cellsInLayers[7]  = 2163;
    m_cellsInLayers[8]  = 2163;
    m_cellsInLayers[9]  = 2194;
    m_cellsInLayers[10] = 2194;
    m_cellsInLayers[11] = 2201;
    m_cellsInLayers[12] = 2201;
    m_cellsInLayers[13] = 2231;
    m_cellsInLayers[14] = 2231;
    m_cellsInLayers[15] = 2238;
    m_cellsInLayers[16] = 2238;
    m_cellsInLayers[17] = 2268;
    m_cellsInLayers[18] = 2268;
    m_cellsInLayers[19] = 2275;
    m_cellsInLayers[20] = 2275;
    m_cellsInLayers[21] = 2312;
    m_cellsInLayers[22] = 2312;
    m_cellsInLayers[23] = 2314;
    m_cellsInLayers[24] = 2314;
    m_cellsInLayers[25] = 2348;
    m_cellsInLayers[26] = 2348;
    m_cellsInLayers[27] = 2356;
    m_cellsInLayers[28] = 2356;
    m_cellsInLayers[29] = 2386;
}

/*****************************************************************/
ProjectionMapProducer::~ProjectionMapProducer()
/*****************************************************************/
{
    m_projectedBins->Write();
    m_unProjected->Write();
    m_output->Close();
}



/*****************************************************************/
void ProjectionMapProducer::initialize()
/*****************************************************************/
{
    m_hgcalNavigator.initialize();
    m_output = TFile::Open("projectionMapping.root", "RECREATE");

    findBinning();

    m_projectedBins = new TH2I("projectedBins", "projectedBins", m_xN, m_xMin, m_xMax, m_yN, m_yMin, m_yMax);
    for(int bx=1;bx<m_projectedBins->GetNbinsX();bx++)
    {
        for(int by=1;by<m_projectedBins->GetNbinsY();by++)
        {
            m_projectedBins->SetBinContent(bx,by,-1);
        }
    }
    m_unProjected = new TH2I("unprojectedBins", "unprojectedBins", m_xN, m_xMin, m_xMax, m_yN, m_yMin, m_yMax);
}



/*****************************************************************/
void ProjectionMapProducer::produce()
/*****************************************************************/
{
    int iz = 1;
    int sec = 1;
    int lay = 1;
    int subsec = 1;
    // fill bin ID in layer 1
    for (int cell=0; cell<m_cellsInLayers[lay-1]; ++cell) 
    {
        const HGCEEDetId id(HGCEE,iz,lay,sec,subsec,cell);
        if (!m_hgcalNavigator.valid(id)) continue;
        const GlobalPoint& point = m_hgcalNavigator.geometry()->getPosition(id);
        double x = point.x();
        double y = point.y();
        int bx = m_projectedBins->GetXaxis()->FindBin(x);
        int by = m_projectedBins->GetYaxis()->FindBin(y);

        m_projectedBins->SetBinContent(bx,by,cell);
    }

    // project all layers onto layer 1
    for (int lay=1; lay<=30; ++lay)
    {
        for (int cell=0; cell<m_cellsInLayers[lay-1]; ++cell) 
        {
            const HGCEEDetId id(HGCEE,iz,lay,sec,subsec,cell);
            if (!m_hgcalNavigator.valid(id)) continue;
            const GlobalPoint& point = m_hgcalNavigator.geometry()->getPosition(id);
            double x = point.x();
            double y = point.y();
            int bx = m_projectedBins->GetXaxis()->FindBin(x);
            int by = m_projectedBins->GetYaxis()->FindBin(y);

            int cellid = m_projectedBins->GetBinContent(bx,by);
            if(cellid!=-1)
            {

                m_cellToBin[make_pair(lay,cell)] = cellid;
                m_binToCell[cellid].push_back(make_pair(lay,cell));
            }
            else
            {
                double drmin = 99999.;
                double xc = m_projectedBins->GetXaxis()->GetBinCenter(bx);
                double yc = m_projectedBins->GetYaxis()->GetBinCenter(by);
                for(int bbx=1;bbx<m_projectedBins->GetNbinsX();bbx++)
                {
                    for(int bby=1;bby<m_projectedBins->GetNbinsY();bby++)
                    {
                        if(m_projectedBins->GetBinContent(bbx,bby)==-1) continue;
                        double xx = m_projectedBins->GetXaxis()->GetBinCenter(bbx);
                        double yy = m_projectedBins->GetYaxis()->GetBinCenter(bby);
                        double dx = (xx-x);
                        double dy = (yy-y);
                        double dr = sqrt(dx*dx+dy*dy);
                        if(dr<drmin)
                        {
                            drmin = dr;
                            cellid = m_projectedBins->GetBinContent(bbx,bby);
                            xc = m_projectedBins->GetXaxis()->GetBinCenter(bbx);
                            yc = m_projectedBins->GetYaxis()->GetBinCenter(bby);
                        }
                    }
                }
                if(point.eta()>1.6 && (fabs(x-xc)>0.475 || fabs(y-yc)>0.475))
                {
                    cout<<cell<<":"<<x-xc<<","<<y-yc<<" (eta="<<point.eta()<<",phi="<<point.phi()<<",layer="<<lay<<")\n";
                    cout<<"   "<<cell<<" -> "<<cellid<<"\n";
                }
                m_unProjected->Fill(x,y);
                m_cellToBin[make_pair(lay,cell)] = cellid;
                m_binToCell[cellid].push_back(make_pair(lay,cell));
                //cout<<"Cannot project cell "<<cell<<" layer "<<lay<<" on layer 1\n";
                //cout<<" Project on "<<cellid<<" instead (dr = "<<drmin<<"\n";
            }
        }
    }

    fstream output("projectionMapping.txt", std::fstream::out);
    if(!output.is_open()) cout<<"Cannot open output file\n";
    for(auto cell_bin : m_cellToBin)
    {
        int layer = cell_bin.first.first;
        int cell  = cell_bin.first.second;
        int id = cell_bin.second;
        output<<layer<<" "<<cell<<" "<<id<<"\n";
    }
    output.close();

    //for(unsigned id=0;id<m_binToCell.size();id++)
    //{
    //    cout<<id<<" -> "<<m_binToCell[id].size()<<"(";
    //    for(unsigned b=0;b<m_binToCell[id].size();b++)
    //    {
    //        cout<<m_binToCell[id][b].first<<","<<m_binToCell[id][b].second<<" ";
    //    }
    //    cout<<")\n";
    //}

}



/*****************************************************************/
void ProjectionMapProducer::findBinning()
/*****************************************************************/
{
    vector< pair<double,double> > cellsXY;

    int iz = 1;
    int sec = 1;
    int lay = 1;
    int subsec = 1;
    // fill cell positions in the first layer
    for (int cell=0; cell<m_cellsInLayers[lay-1]; ++cell) 
    {
        const HGCEEDetId id(HGCEE,iz,lay,sec,subsec,cell);
        if (!m_hgcalNavigator.valid(id)) continue;
        const GlobalPoint& point = m_hgcalNavigator.geometry()->getPosition(id);
        double x = point.x();
        double y = point.y();

        cellsXY.push_back(make_pair(x,y));
    }

    // find min, max, and bin size
    double drMin = 99999.;
    double xMin  = 999999.;
    double xMax  = 0.;
    double yMin  = 999999.;
    double yMax  = -999999.;
    for(unsigned i=0;i<cellsXY.size();i++)
    {
        double xi = cellsXY[i].first;
        double yi = cellsXY[i].second;
        if(xi<xMin) xMin = xi;
        if(xi>xMax) xMax = xi;
        if(yi<yMin) yMin = yi;
        if(yi>yMax) yMax = yi;

        for(unsigned j=i+1;j<cellsXY.size();j++)
        {
            double xj = cellsXY[j].first;
            double yj = cellsXY[j].second;
            double dx = fabs(xj-xi);
            double dy = fabs(yj-yi);
            double dr = sqrt(dx*dx+dy*dy);
            if( dr<drMin && dr>0. ) 
            {
                drMin = dr;
            }
        }
    }
    xMin -= drMin/2.;
    xMax += drMin/2.;
    yMin -= drMin/2.;
    yMax += drMin/2.;
    int xN = (xMax-xMin)/drMin;
    int yN = (yMax-yMin)/drMin;
    xMin -= drMin*(int)(xN/10);
    xMax += drMin*(int)(xN/10);
    yMin -= drMin*(int)(yN/10);
    yMax += drMin*(int)(yN/10);

    cout<<(xMax-xMin)/drMin<<"\n";
    m_xN = (xMax-xMin)/drMin;
    m_xMin = xMin;
    m_xMax = xMax;

    m_yN = (yMax-yMin)/drMin;
    m_yMin = yMin;
    m_yMax = yMax;

    cout<<"x min="<<xMin<<", max="<<xMax<<"\n";
    cout<<"y min="<<yMin<<", max="<<yMax<<"\n";
    cout<<"bin size = "<<drMin<<"\n";

}
