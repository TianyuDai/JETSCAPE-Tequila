#include <iostream>
#include <fstream>
std::ifstream fin("../build/pp2760_colored_brick/sigmaGen"); 
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
// #include "JetScapePartonReader.h"
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
	const int nEvents = 10000; 

	vector <double> pTBin{64., 74., 84., 97., 114., 133., 153., 174., 196., 220, 245, 272, 300};
    // vector <double> pTHatBin{5., 10., 50., 100., 200., 400., 600., 800., 1000., 1380.};
    vector <double> pTHatBin{20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.};
	// vector <double> pTHatBin{5., 10., 20., 40., 60., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
	// vector <double> pTHatBin{160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.}; 
    vector <double> sigmaGen(pTHatBin.size()-1);  
    vector <double> sigmaErr(pTHatBin.size()-1);  
	// std::vector <double> jet_cs(pTBin.size()-1, 0.), err(pTBin.size()-1, 0.), sqSum(pTBin.size()-1, 0.); 
    //std::vector <int> jet_ct(pTBin.size()-1);
    for (unsigned int iSigma = 0; iSigma < pTHatBin.size() - 1; iSigma++)
    {
        double pTHat; 
        fin >> pTHat; 
        if (pTHat == pTHatBin[iSigma])
            fin >> sigmaGen[iSigma] >> sigmaErr[iSigma]; 
        else std::cout << pTHat << " pTHatBin boundary does not match! \n"; 
    } 

    double jetpTMin = 50., jetRadius = 0.2, jetAbsRapMax = 2.; 
    double deltaRap = jetAbsRapMax * 2.; 
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
    vector <fjcore::PseudoJet> fjInputs; 
    fjcore::Selector select_rapidity = fjcore::SelectorAbsRapMax(jetAbsRapMax); 

	std::ofstream jet_output ("/home/td115/research/Result/JETSCAPE2.0/general/pp/pp2760_colored_hadron_brickHydro_inclusiveJet_R0.2_pTMin50.txt");
	// std::ofstream jet_output ("test.txt");
   	std::vector <double> jet_cs(pTBin.size()-1, 0.), jet_cs_err(pTBin.size()-1, 0.); 
    std::vector<std::vector <int>> jet_ct(pTHatBin.size()-1);
    std::vector<std::vector <int>> jet_ct_sq(pTHatBin.size()-1);
	for (unsigned int iBin = 0; iBin < pTHatBin.size() - 1; iBin++)
	{
		auto reader=make_shared<JetScapeReaderAscii>(("../build/pp2760_colored_brick/"+std::to_string(pTHatBin[iBin])+".dat").c_str());  
		jet_ct[iBin] = std::vector<int>(pTBin.size()-1, 0); 
		jet_ct_sq[iBin] = std::vector<int>(pTBin.size()-1, 0);
                // pT_sq[iBin] = std::vector<int>(pTBin.size()-1, 0); 
		while (!reader->Finished())
		{
		        std::vector<int> jet_ct_s(pTBin.size()-1, 0); 
			reader->Next();
			fjInputs.resize(0);          
			vector <fjcore::PseudoJet> inclusiveJets, sortedJets;
			// fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef); 
			fjcore::ClusterSequence hadroClustSeq(reader->GetHadronsForFastJet(), jetDef); 
			inclusiveJets = hadroClustSeq.inclusive_jets(jetpTMin); 
			// std::cout << inclusiveJets.size() << endl; 
			vector <fjcore::PseudoJet> selected_jets = select_rapidity(inclusiveJets); 
			sortedJets = fjcore::sorted_by_pt(selected_jets); 
			for (unsigned int iJet = 0; iJet < sortedJets.size(); iJet++)
			{
                            // if (iJet > 0) break; 
				for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
					if (sortedJets[iJet].perp() > pTBin[ipT] && sortedJets[iJet].perp() < pTBin[ipT+1])
				    {
						// std::cout << "a jet" << endl; 
				    	jet_ct_s[ipT]++; 
				        break; 
				    }
			}
                        for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
                        {
                            jet_ct[iBin][ipT] += jet_ct_s[ipT]; 
                            jet_ct_sq[iBin][ipT] += jet_ct_s[ipT] * jet_ct_s[ipT]; 
                        }
                        // jet_err[iBin][ipT]
		}
		// double sigma = sigmaGen[iBin] / nEvents; 
		for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
		{
			jet_cs[ipT] += double(jet_ct[iBin][ipT]) * sigmaGen[iBin] / nEvents; 
                        // jet_cs_sq[ipT] += jet_ct_sq[iBin][ipT] * pow(sigmaGen[iBin], 2) / nEvents; 
                        double err_ipT; 
                        err_ipT = sqrt((double(jet_ct_sq[iBin][ipT])/nEvents - pow(double(jet_ct[iBin][ipT])/nEvents, 2))/nEvents);
                        // std::cout << double(jet_ct_sq[iBin][ipT])/nEvents - pow(double(jet_ct[iBin][ipT])/nEvents, 2) << "\n"; 
                        jet_cs_err[ipT] += sqrt(sigmaErr[iBin]*sigmaErr[iBin] * (double(jet_ct[iBin][ipT])/nEvents)*(double(jet_ct[iBin][ipT])/nEvents) + sigmaGen[iBin]*sigmaGen[iBin]*err_ipT*err_ipT); 
                        // std::cout << iBin << " " << ipT << " " << err_ipT << " " << jet_cs_err[ipT] << "\n"; 
			// sqSum[ipT] += jet_ct[iBin][ipT] * pow(sigmaGen[iBin]/nEvents, 2);
		}
		reader->Close(); 
	}
	for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
    {
        		// err[ipT] = jet_cs[ipT] / sqrt(pow(jet_cs[ipT], 2) / sqSum[ipT]); 
			jet_output << (pTBin[ipT] + pTBin[ipT+1]) / 2 << " " << jet_cs[ipT] / (pTBin[ipT+1] - pTBin[ipT]) / 4 * 1.e6 << " " << jet_cs_err[ipT] / (pTBin[ipT+1] - pTBin[ipT]) / 4 * 1.e6  << endl; 
	}

	cout<<"Finished!"<<endl; 
}
