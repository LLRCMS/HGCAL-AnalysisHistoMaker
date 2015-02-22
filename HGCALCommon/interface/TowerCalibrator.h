
#ifndef TOWERCALIBRATOR_H
#define TOWERCALIBRATOR_H

#include "AnHiMaHGCAL/HGCALCommon/interface/Tower.h"

#include <vector>


namespace AnHiMa
{

    class TowerCalibrator
    {
        public:
            TowerCalibrator();
            ~TowerCalibrator();

            void calibrate(Tower& tower, bool halfLayers=false);
            double calibratedEnergy(double energy, double eta, int layer, int subdet=3);

        private:
            // EE
            double m_mipValueInGeV_EE;
            double m_coeff_a_EE;
            double m_coeff_b_EE;
            double m_coeff_c_EE;
            std::vector<double> m_weights_EE;
            // HEF
            double m_mipValueInGeV_HEF;
            double m_coeff_a_HEF;
            double m_coeff_b_HEF;
            double m_coeff_c_HEF;
            std::vector<double> m_weights_HEF;
            // HEB
            double m_mipValueInGeV_HEB;
            double m_coeff_a_HEB;
            double m_coeff_b_HEB;
            double m_coeff_c_HEB;
            std::vector<double> m_weights_HEB;

            // 1/2 layer calibration
            std::vector<double> m_halfLayersCalib;
    };

};

#endif
