
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

            void calibrate(Tower& tower);
            double calibratedEnergy(double energy, double eta, int layer);

        private:
            double m_mipValueInGeV;
            double m_coeff_a;
            double m_coeff_b;
            double m_coeff_c;
            std::vector<double> m_weights;

    };

};

#endif
