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
ofstream fout("../Nersc_submit/tasks_list_AA_event_cut20.txt"); 

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
    vector <int> iTask{123, 151, 15, 187, 193, 229, 232, 245, 260, 282, 293, 310, 329, 349, 364, 380, 382, 391, 396, 413, 426, 428, 437, 43, 457, 458, 478, 480, 487, 521, 533, 549, 556, 572, 577, 579, 587, 595, 596, 596, 607, 609, 616, 633, 672, 675, 67, 697, 698, 711, 736, 744, 745, 754, 758, 761, 762, 765, 766, 775, 785, 787, 78, 791, 794, 799, 818, 819, 827, 829, 831, 840, 876, 877, 878, 87, 884, 889, 897, 912, 937}; 

    fout << "#!/usr/bin/env bash\n\n"; 
    for (size_t iBin = 0; iBin < pTHatBin.size()-1; iBin++)
    {
        double pTHatMin = pTHatBin[iBin]; 
        double pTHatMax = pTHatBin[iBin+1]; 
        for (int i = 0; i < 100; i++)
            fout << "cd /global/homes/t/td115/JETSCAPE3.0-Tequila/build && python3 AA_wrapper.py --task "+std::to_string(i)+" --pTHatMin "+std::to_string(pTHatMin)+" --pTHatMax "+std::to_string(pTHatMax) << "\n";  
    }/*
    for (size_t iBin = 0; iBin < pTHatBin.size()-1; iBin++)
    {
        double pTHatMin = pTHatBin[iBin]; 
        double pTHatMax = pTHatBin[iBin+1]; 
        for (int i = 0; i < 100; i++)
            fout << "cd /global/homes/t/td115/JETSCAPE2.0-AAevent/JETSCAPE2.0-Tequila/build && ./pp_charged_hadron_cross_section --task "+std::to_string(i)+" --pTHatMin "+std::to_string(pTHatMin)+" --pTHatMax "+std::to_string(pTHatMax) << "\n";  
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
