// compared to CMS paper - arXiv:1202.2554
#include <iostream>
#include <fstream>
std::ifstream fin("../build/pp2760_colorless_brick/sigmaGen"); 
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
    double total_cross_section = 78.1179; 
    double d_eta = 2., d_phi = 2.*M_PI; 
    vector <double> pTBin{7.2, 10.8, 14.4, 21.6, 28.8, 38.4, 48., 67.2, 86.4, 112.2, 120.};
    vector <double> pTHatBin{10., 20., 50., 80., 110., 160., 210., 260., 310., 400., 500., 600., 800., 1000., 1380.};
    vector <double> sigmaGen(pTHatBin.size()-1);  
    vector <double> sigmaErr(pTHatBin.size()-1);  
    for (unsigned int iSigma = 0; iSigma < pTHatBin.size() - 1; iSigma++)
    {
        double pTHat; 
        fin >> pTHat; 
        if (pTHat == pTHatBin[iSigma])
            fin >> sigmaGen[iSigma] >> sigmaErr[iSigma]; 
        else std::cout << pTHat << " pTHatBin boundary does not match! \n"; 
    }


    std::ofstream jet_output ("/home/td115/research/Result/JETSCAPE2.0/general/AA/PbPb2760_colorless_hadron_brickHydro_chargedHadron`.txt");

    std::vector <double> hadron_cs(pTBin.size()-1, 0.), hadron_cs_err(pTBin.size()-1, 0.), pTAvg(pTBin.size()-1, 0.), total_sigma(pTBin.size()-1, 0.); 
    std::vector<std::vector <int>> hadron_ct(pTHatBin.size()-1);
    std::vector<std::vector <int>> hadron_ct_sq(pTHatBin.size()-1);
    for (unsigned int iBin = 0; iBin < pTHatBin.size() - 1; iBin++)
    {
        auto reader=make_shared<JetScapeReaderAscii>(("../build/pp2760_colorless_brick/"+std::to_string(pTHatBin[iBin])+".dat").c_str());  
	hadron_ct[iBin] = std::vector<int>(pTBin.size()-1, 0); 
	hadron_ct_sq[iBin] = std::vector<int>(pTBin.size()-1, 0);
	std::vector<double> pTSum(pTBin.size()-1, 0.);
	while (!reader->Finished())
	{
            std::vector<int> hadron_ct_s(pTBin.size()-1, 0); 
            reader->Next();
            auto hadrons = reader->GetHadrons(); 
            auto partons = reader->GetPartons();
	    auto pdghelper = JetScapeParticleBase::InternalHelperPythia.particleData; 
            for (unsigned int iHadron = 0; iHadron < hadrons.size(); iHadron++)
            {
                if (hadrons[iHadron]->pt() < pTBin[0]) continue; 
                if (fabs(hadrons[iHadron]->eta()) > 1.) continue; 
                auto ID = hadrons[iHadron]->pid(); 
                auto charge = pdghelper.charge( ID );
                if (charge == 0) continue; 
		for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
		    if (hadrons[iHadron]->pt() > pTBin[ipT] && hadrons[iHadron]->pt() <= pTBin[ipT+1])
		    {
		        hadron_ct_s[ipT]++; 
                        pTSum[ipT] += hadrons[iHadron]->pt(); 
		        break; 
		    }
            }
            for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
            {
                hadron_ct[iBin][ipT] += hadron_ct_s[ipT]; 
                hadron_ct_sq[iBin][ipT] += hadron_ct_s[ipT] * hadron_ct_s[ipT]; 
            }
	}
	for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
	{
            if (hadron_ct[iBin][ipT] > 0)
            {
               pTSum[ipT] /= hadron_ct[iBin][ipT]; 
               total_sigma[ipT] += sigmaGen[iBin]; 
            }
            pTAvg[ipT] += pTSum[ipT] * sigmaGen[iBin]; 
	    hadron_cs[ipT] += double(hadron_ct[iBin][ipT]) / nEvents * sigmaGen[iBin];
            std::cout << "iBin " << iBin << " ipT " << ipT << " cross section " << double(hadron_ct[iBin][ipT]) / nEvents * sigmaGen[iBin] << "\n";  
            double err_ipT; 
            err_ipT = sqrt((double(hadron_ct_sq[iBin][ipT])/nEvents - pow(double(hadron_ct[iBin][ipT])/nEvents, 2))/nEvents);
            hadron_cs_err[ipT] += sqrt(sigmaErr[iBin]*sigmaErr[iBin] * (double(hadron_ct[iBin][ipT])/nEvents)*(double(hadron_ct[iBin][ipT])/nEvents) + sigmaGen[iBin]*sigmaGen[iBin]*err_ipT*err_ipT); 
	}
	reader->Close(); 
    }
    for (unsigned int ipT = 0; ipT < pTBin.size()-1; ipT++)
	jet_output << pTAvg[ipT]/total_sigma[ipT] << " " << hadron_cs[ipT] / (pTBin[ipT+1] - pTBin[ipT]) / d_eta / d_phi / pTAvg[ipT] * total_sigma[ipT] / total_cross_section << " " << hadron_cs_err[ipT] / (pTBin[ipT+1] - pTBin[ipT]) / d_eta /d_phi / pTAvg[ipT] * total_sigma[ipT] / total_cross_section << endl; 

    cout<<"Finished!"<<endl; 
}
