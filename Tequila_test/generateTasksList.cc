#include <fstream>
#include <iostream>
#include <time.h>
// JetScape Framework includes ...
#include "JetScape.h"
#include "JetEnergyLoss.h"
#include "JetEnergyLossManager.h"
#include "JetScapeWriterStream.h"
#ifdef USE_HEPMC
#include "JetScapeWriterHepMC.h"
#endif


// User modules derived from jetscape framework clasess
#include "TrentoInitial.h"
#include "Tequila.h"
#include "Martini.h"
#include "Brick.h"
#include "PGun.h"
#include "PythiaGun.h"
#include "PartonPrinter.h"
#include "HydroFromFile.h"

#include "tinyxml2.h"
#include <chrono>
#include <thread>

#define hbarc 0.197327053

using namespace std;
ofstream fout("tasks_list_AA_event.txt"); 

using namespace Jetscape;

// Forward declaration
void Show();

// -------------------------------------

int main(int argc, char** argv)
{
    clock_t t; t = clock();
    time_t start, end; time(&start);

    JetScapeLogger::Instance()->SetInfo(true);
    JetScapeLogger::Instance()->SetDebug(false);
    JetScapeLogger::Instance()->SetRemark(false);
    JetScapeLogger::Instance()->SetVerboseLevel(0);
  
    Show();

    cout<<endl;
   
    vector <double> pTHatBin{5., 10., 20., 40., 60., 80., 100.};
    vector <double> cut_list{1.}; 
    vector <double> coef_list{1.}; 

    fout << "#!/usr/bin/env bash\n\n"; 
    for (size_t iBin = 0; iBin < pTHatBin.size()-1; iBin++)
    {
        double pTHatMin = pTHatBin[iBin]; 
        double pTHatMax = pTHatBin[iBin+1]; 
      for (size_t iCut = 0; iCut < cut_list.size(); iCut++)
      {
          double cut = cut_list[iCut]; 
          double coef = coef_list[iCut]; 
          for (int i = 0; i < 10; i++)
              fout << "cd /home/td115/research/JETSCAPE3.0/JETSCAPE3.0-Tequila/build && python3 AA_wrapper.py --task "+std::to_string(i)+" --pTHatMin "+std::to_string(pTHatMin)+" --pTHatMax "+std::to_string(pTHatMax)+" --cut "+std::to_string(cut)+" --coef "+std::to_string(coef) << "\n";  
      }
    }
/*
    for (size_t iBin = 0; iBin < pTHatBin.size()-1; iBin++)
    {
        double pTHatMin = pTHatBin[iBin]; 
        double pTHatMax = pTHatBin[iBin+1]; 
        for (int i = 0; i < 100; i++)
            fout << "cd /global/homes/t/td115/JETSCAPE3.0-Tequila/build && ./chargedHadronCrossSection --task "+std::to_string(i)+" --pTHatMin "+std::to_string(pTHatMin)+" --pTHatMax "+std::to_string(pTHatMax) << "\n";  
    }
  
 */
  INFO_NICE<<"Finished!"; 
  cout<<endl;

  t = clock() - t;
  time(&end);
  printf ("CPU time: %f seconds.\n",((float)t)/CLOCKS_PER_SEC);
  printf ("Real time: %f seconds.\n",difftime(end,start));

  return 0;
}

// -------------------------------------

void Show()
{
  INFO_NICE<<"------------------------------------";
  INFO_NICE<<"| Brick Test JetScape Framework ... |";
  INFO_NICE<<"------------------------------------";
  INFO_NICE;
}
