

#include <sstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/ProjectionMapping/interface/ProjectionMapProducer.h"

using namespace std;
using namespace AnHiMa;


const int ProjectionMapProducer::NCELLS = 2145;

/*****************************************************************/
ProjectionMapProducer::ProjectionMapProducer():
    m_output(0),
    m_projectedBins(0),
    m_unProjected(0),
    m_binToCell(NCELLS),
    m_xN(0),
    m_yN(0),
    m_xMin(0.),
    m_xMax(0.),
    m_yMin(0.),
    m_yMax(0.)
/*****************************************************************/
{
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
    for (int cell=0; cell<NCELLS; ++cell) 
    {
        const HGCEEDetId id(HGCEE,iz,lay,sec,subsec,cell);
        if (!m_hgcalNavigator.valid(id)) continue;
        const GlobalPoint& point = m_hgcalNavigator.geometry()->getPosition(id);
        double x = point.x();
        double y = point.y();
        int bx = m_projectedBins->GetXaxis()->FindBin(x);
        int by = m_projectedBins->GetYaxis()->FindBin(y);
        //cout<<cell<<":"<<x<<","<<y<<" -> "<<bx<<","<<by<<"\n";
        m_projectedBins->SetBinContent(bx,by,cell);
    }

    // project all layers onto layer 1
    for (int lay=1; lay<=30; ++lay)
    {
        for (int cell=0; cell<NCELLS; ++cell) 
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
                        }
                    }
                }
                m_unProjected->Fill(x,y);
                cout<<"Cannot project cell "<<cell<<" layer "<<lay<<" on layer 1\n";
                cout<<" Project on "<<cellid<<" instead (dr = "<<drmin<<"\n";
            }
        }
    }

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
    // fill cell positions
    for (int cell=0; cell<NCELLS; ++cell) 
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
