/**
 *  @file  IObjectFactory.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    07/29/2014
 *
 *  @internal
 *     Created :  07/29/2014
 * Last update :  07/29/2014 02:55:59 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */



#ifndef IOBJECTFACTORY_H
#define IOBJECTFACTORY_H

#include "IEvent.h"

#include "TChain.h"

namespace AnHiMa
{
    template<class T> class IObjectFactory
    {
        public:
            IObjectFactory() {};
            virtual ~IObjectFactory() {};

            const T& data() const {return m_data;}
            virtual void initialize(IEvent*, TChain*) = 0;


        protected:
            virtual void connectVariables(TChain*) = 0;
            virtual void update() = 0;

            T m_data;

    };

}


#endif
