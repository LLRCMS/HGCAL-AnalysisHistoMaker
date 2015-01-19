
#ifndef TPARTICLEPAIR_H
#define TPARTICLEPAIR_H

#include <vector>
#include "TLorentzVector.h"


namespace AnHiMa
{
    template<class T, class Tbase=TLorentzVector> class TParticlePair : public Tbase
    {
        public:
            TParticlePair(const TParticlePair<T,Tbase>& p):Tbase(p)
            {
                m_constituents.push_back(&p[0]);
                m_constituents.push_back(&p[1]);
            }
            TParticlePair(const T& p1, const T& p2):Tbase(p1)
            {
                (*this) += static_cast<const Tbase&>(p2); 
                m_constituents.push_back(&p1);
                m_constituents.push_back(&p2);
            }
            ~TParticlePair() {};


            const T& operator[](std::size_t i) const {return *(m_constituents[i]);}
            const T& constituent(std::size_t i) const {return *(m_constituents[i]);}

        private:
            std::vector<const T*> m_constituents;
    };
}

#endif
