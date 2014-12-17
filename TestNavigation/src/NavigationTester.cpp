

#include <sstream>

#include "TVector2.h"

#include "AnHiMaHGCAL/TestNavigation/interface/NavigationTester.h"

using namespace std;
using namespace AnHiMa;



/*****************************************************************/
NavigationTester::NavigationTester():
    m_output(0),
    m_cellsSource(0),
    m_cells_1_0(0),
    m_cells_2_0(0),
    m_cells_3_0(0),
    m_cells_m1_0(0),
    m_cells_m2_0(0),
    m_cells_m3_0(0)
/*****************************************************************/
{
}

/*****************************************************************/
NavigationTester::~NavigationTester()
/*****************************************************************/
{
    m_cellsSource->Write();
    m_cells_1_0->Write();
    m_cells_2_0->Write();
    m_cells_3_0->Write();
    m_cells_m1_0->Write();
    m_cells_m2_0->Write();
    m_cells_m3_0->Write();

    m_output->Close();
}



/*****************************************************************/
void NavigationTester::initialize()
/*****************************************************************/
{
    m_hgcalNavigator.initialize();
    m_output = TFile::Open("navigationTest.root", "RECREATE");

}



/*****************************************************************/
void NavigationTester::leftRightTest()
/*****************************************************************/
{
    vector<double> xs;
    vector<double> xs_1_0;
    vector<double> xs_2_0;
    vector<double> xs_3_0;
    vector<double> xs_m1_0;
    vector<double> xs_m2_0;
    vector<double> xs_m3_0;

    vector<double> ys;
    vector<double> ys_1_0;
    vector<double> ys_2_0;
    vector<double> ys_3_0;
    vector<double> ys_m1_0;
    vector<double> ys_m2_0;
    vector<double> ys_m3_0;

    int iz = 1;
    int sec = 1;
    int lay = 1;
    int subsec = 1;
    //for (int sub=0; sub<=1; ++sub) 
    //{
    //    int subsec = (sub==0 ? -1 : 1);
        for (int cell=0; cell<700; cell+=5) 
        //for (int cell=0; cell<700; cell+=1)
        {
            const HGCEEDetId id(HGCEE,iz,lay,sec,subsec,cell);
            if (!m_hgcalNavigator.valid(id)) continue;
            std::cout<<"toto\n";
            const GlobalPoint& point = m_hgcalNavigator.geometry()->getPosition(id);
            double x = point.x();
            double y = point.y();
            //double z = point.z();
            //double eta = point.eta();
            //double phi = point.phi();
            //if(fabs(phi)>0.175) continue;
            std::cout<<"1\n";

            // navigate
            DetId id_1_0   = m_hgcalNavigator.topology()->offsetByFix(id, 1,  0);
            DetId id_2_0   = m_hgcalNavigator.topology()->offsetByFix(id, 2,  0);
            DetId id_3_0   = m_hgcalNavigator.topology()->offsetByFix(id, 3,  0);
            DetId id_m1_0  = m_hgcalNavigator.topology()->offsetByFix(id, -1,  0);
            DetId id_m2_0  = m_hgcalNavigator.topology()->offsetByFix(id, -2,  0);
            DetId id_m3_0  = m_hgcalNavigator.topology()->offsetByFix(id, -3,  0);

            const GlobalPoint& point_1_0  = m_hgcalNavigator.geometry()->getPosition(id_1_0);
            const GlobalPoint& point_2_0  = m_hgcalNavigator.geometry()->getPosition(id_2_0);
            const GlobalPoint& point_3_0  = m_hgcalNavigator.geometry()->getPosition(id_3_0);
            const GlobalPoint& point_m1_0 = m_hgcalNavigator.geometry()->getPosition(id_m1_0);
            const GlobalPoint& point_m2_0 = m_hgcalNavigator.geometry()->getPosition(id_m2_0);
            const GlobalPoint& point_m3_0 = m_hgcalNavigator.geometry()->getPosition(id_m3_0);

            double x_1_0  = point_1_0.x();
            double x_2_0  = point_2_0.x();
            double x_3_0  = point_3_0.x();
            double x_m1_0 = point_m1_0.x();
            double x_m2_0 = point_m2_0.x();
            double x_m3_0 = point_m3_0.x();

            double y_1_0  = point_1_0.y();
            double y_2_0  = point_2_0.y();
            double y_3_0  = point_3_0.y();
            double y_m1_0 = point_m1_0.y();
            double y_m2_0 = point_m2_0.y();
            double y_m3_0 = point_m3_0.y();

            xs     .push_back(x);
            xs_1_0 .push_back(x_1_0);
            xs_2_0 .push_back(x_2_0);
            xs_3_0 .push_back(x_3_0);
            xs_m1_0.push_back(x_m1_0);
            xs_m2_0.push_back(x_m2_0);
            xs_m3_0.push_back(x_m3_0);

            ys     .push_back(y);
            ys_1_0 .push_back(y_1_0);
            ys_2_0 .push_back(y_2_0);
            ys_3_0 .push_back(y_3_0);
            ys_m1_0.push_back(y_m1_0);
            ys_m2_0.push_back(y_m2_0);
            ys_m3_0.push_back(y_m3_0);

        }
    //}

    m_cellsSource = new TGraph(xs.size(), &xs[0], &ys[0]);
    m_cells_1_0   = new TGraph(xs_1_0.size(), &xs_1_0[0], &ys_1_0[0]);
    m_cells_2_0   = new TGraph(xs_2_0.size(), &xs_2_0[0], &ys_2_0[0]);
    m_cells_3_0   = new TGraph(xs_3_0.size(), &xs_3_0[0], &ys_3_0[0]);
    m_cells_m1_0  = new TGraph(xs_m1_0.size(), &xs_m1_0[0], &ys_m1_0[0]);
    m_cells_m2_0  = new TGraph(xs_m2_0.size(), &xs_m2_0[0], &ys_m2_0[0]);
    m_cells_m3_0  = new TGraph(xs_m3_0.size(), &xs_m3_0[0], &ys_m3_0[0]);

    m_cellsSource->SetName("cellsSource");
    m_cells_1_0->SetName("cells_1_0");
    m_cells_2_0->SetName("cells_2_0");
    m_cells_3_0->SetName("cells_3_0");
    m_cells_m1_0->SetName("cells_m1_0");
    m_cells_m2_0->SetName("cells_m2_0");
    m_cells_m3_0->SetName("cells_m3_0");

}



