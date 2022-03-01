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
#include "JetScapeReader.h"
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


  double jetpTMin = 5., jetRadius = 0.6, jetAbsRapMax = 1., partonpTMin = 2.;

  fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
  vector <fjcore::PseudoJet> fjInputs;
  fjcore::Selector select_rapidity = fjcore::SelectorAbsRapMax(jetAbsRapMax);
  fjcore::Selector select_pt = fjcore::SelectorPtMin(partonpTMin); 
  fjcore::Selector select_both = select_rapidity && select_pt; 

  while (!reader->Finished())
  {
    reader->Next();

    auto hadrons_for_jet = reader->GetHadronsForFastJet();
    // auto hadrons = reader->GetHadrons();
    // for (unsigned int i = 0; i < )
    

    if(hadrons_for_jet.size() > 0)
    {
      ++SN;
      dist_output << "#" << "\t"
            << "E"  <<"\t"
            << "Px"  << "\t"
            << "Py"  << "\t"
            << "Pz"  << "\t"
            << "Pt" <<  "\t"
            << "eta" << "\t" << "Phi" << "\n";

      fjInputs.resize(0);
      vector <fjcore::PseudoJet> inclusiveJets, sortedJets;


      fjcore::ClusterSequence clustSeq(hadrons_for_jet, jetDef);
      inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
      vector <fjcore::PseudoJet> selected_jets = select_both(inclusiveJets);
      sortedJets = fjcore::sorted_by_pt(selected_jets);
      
      for (unsigned int iJet = 0; iJet < sortedJets.size(); iJet++)
      {
          dist_output << iJet
              << particleSeperator << sortedJets[iJet].e()
              << particleSeperator << sortedJets[iJet].px()
              << particleSeperator << sortedJets[iJet].py()
              << particleSeperator << sortedJets[iJet].pz()
              << particleSeperator << sortedJets[iJet].perp() 
              << particleSeperator << sortedJets[iJet].eta() 
              << particleSeperator << sortedJets[iJet].phi() <<"\n"; 
          // std::cout << "jet " << iJet << " " << sortedJets[iJet].perp() << "\n"; }
    }
   }
  }
  // Write the final cross section and error if requested by using header v2
  if (writeHeaderV2) {
    // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
    dist_output << "#"
        << "\t" << "sigmaGen\t" << reader->GetSigmaGen()
        << "\t" << "sigmaErr\t" << reader->GetSigmaErr()
        << "\n";
  }
  reader->Close();
}
