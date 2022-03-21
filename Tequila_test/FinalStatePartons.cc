/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapePartonReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
using namespace fjcore;

using namespace Jetscape;

// You could overload here and then simply use ofstream << p;
// ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);


// -------------------------------------

int main(int argc, char** argv)
{

  // JetScapeLogger::Instance()->SetInfo(false);
  JetScapeLogger::Instance()->SetDebug(false);
  JetScapeLogger::Instance()->SetRemark(false);
  // //SetVerboseLevel (9 a lot of additional debug output ...)
  // //If you want to suppress it: use SetVerboseLevle(0) or max  SetVerboseLevel(9) or 10
  JetScapeLogger::Instance()->SetVerboseLevel(0);

  // Whether to write the new header (ie. v2), including xsec info.
  // To enable, pass anything as the third argument to enable this option.
  // Default: disabled.
  bool writeHeaderV2 = false;
  if (argc > 3) {
    writeHeaderV2 = static_cast<bool>(argv[3]);
    std::cout << "NOTE: Writing header v2, and final cross section and error at EOF.\n";
  }

  // The seperator between particles depends on the header.
  std::string particleSeperator = " ";
  if (!writeHeaderV2) {
    particleSeperator = "\t";
  }

  auto reader=make_shared<JetScapeReaderAscii>(argv[1]);
  std::ofstream dist_output (argv[2]); //Format is SN, PID, E, Px, Py, Pz, Eta, Phi
  int SN=0, TotalPartons=0;

  auto pdghelper = JetScapeParticleBase::InternalHelperPythia.particleData;

  while (!reader->Finished())
  {
    reader->Next();

    // cout<<"Analyze current event: "<<reader->GetCurrentEvent()<<endl;

    //dist_output<<"Event "<< reader->GetCurrentEvent()+1<<endl;
    auto partons = reader->GetPartons();
    cout<<"Number of partons is: " << partons.size() << endl;

    if (partons.size() > 0)
    {
      ++SN;
        dist_output << "#" << "\t"
            << "pid" << "\t"
            << "pstat" << "\t"
            << "e"   << "\t"
            << "px"   << "\t"
            << "Py"  << "\t"
            << "Pz"  << "\t"
            << "Eta" <<  "\t"<< "Phi" << "\n";

      for (unsigned int i=0; i<partons.size(); i++)
      {
        auto ID = partons[i].get()->pid(); 
        auto charge = pdghelper.charge(ID); 
        auto eta = partons[i].get()->eta(); 
        
        dist_output << i
            << particleSeperator << partons[i].get()->pid()
            << particleSeperator << partons[i].get()->pstat()
            << particleSeperator << partons[i].get()->e()
            << particleSeperator << partons[i].get()->px()
            << particleSeperator << partons[i].get()->py()
            << particleSeperator << partons[i].get()->pz()
            << particleSeperator << partons[i].get()->eta()
            << particleSeperator << partons[i].get()->phi();
        
        // v2 drops eta and phi, so only include it for v1
        
        // Finish up
        dist_output << "\n";
      }
    }
  }

  reader->Close();
}
