#ifndef TEQUILAMUTEX_H
#define TEQUILAMUTEX_H

#include <vector>
#include <string>
#include<iostream>
#include <memory>

#include "JetScapeTask.h"
#include "JetScapeModuleMutex.h"

using namespace Jetscape;
using std::shared_ptr;


class TequilaMutex : public JetScapeModuleMutex
{
  public:
    TequilaMutex();
    ~TequilaMutex();
    bool CheckMutex(vector<shared_ptr<JetScapeTask>> modules);  



};


#endif
