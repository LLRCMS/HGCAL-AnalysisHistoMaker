/**
 *  @file  IEvent.h
 *  @brief  
 *
 *
 *  @author  Jean-Baptiste Sauvan <sauvan@llr.in2p3.fr>
 *
 *  @date    12/17/2012
 *
 *  @internal
 *     Created :  12/17/2012
 * Last update :  12/17/2012 05:27:28 PM
 *          by :  JB Sauvan
 *
 * =====================================================================================
 */

#ifndef IEVENT_H
#define IEVENT_H

#include <map>
#include <string>
#include <iostream>
#include <algorithm>

class TChain;
class TBranch;

namespace AnHiMa
{
    class IEvent
    {
        private: 
            typedef void (*CallBack)(void*);

        public:
            IEvent();
            ~IEvent();

            virtual bool passSelection(int selection=0) = 0;
            void connectVariables(TChain* inputChain);
            void registerCallback(void* objectPtr, CallBack function);
            virtual void update();

        //protected:
        //    TChain* inputChain() {return m_inputChain;}

        private:
            TChain* m_inputChain;
            std::vector< std::pair<void*, CallBack> > m_listOfCallBacks;



    };
}

#endif
