#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
// #include "JetScapePartonReader.h"
#include "JetScapeBanner.h"
#include <complex>
// #include "fjcore.hh"

#include <GTL/dfs.h>

using namespace std;
// using namespace fjcore;
using namespace Jetscape;

// You could overload here and then simply use ofstream << p;
// ostream & operator<<(ostream & ostr, const fjcore::PseudoJet & jet);


// -------------------------------------

char* getCmdOptions(char ** begin, char ** end, const std::string & option)
{
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr; 
	return 0; 
}

int main(int argc, char** argv)
{
    char * task = getCmdOptions(argv, argv+argc, "--task"); 
    std::string task_str = ""; 
    int i_task = 0; 
    if (task != NULL) {task_str = task; i_task = std::stoi(task_str);} 

    char * pT_hat_min = getCmdOptions(argv, argv+argc, "--pTHatMin"); 
    std::string pT_hat_min_str = ""; 
    double i_pT_hat_min = 0.; 
    if (pT_hat_min != NULL) {pT_hat_min_str = pT_hat_min; i_pT_hat_min = std::stod(pT_hat_min_str);} 

    char * pT_hat_max = getCmdOptions(argv, argv+argc, "--pTHatMax"); 
    std::string pT_hat_max_str = ""; 
    double i_pT_hat_max = 100.; 
    if (pT_hat_max != NULL) {pT_hat_max_str = pT_hat_max; i_pT_hat_max = std::stod(pT_hat_max_str);} 

    vector <double> pTBin{10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48., 50., 52., 54., 56., 58., 60.};
//    std::complex<double> i(0., 1.); 

            std::vector <double> hadron_pt(pTBin.size()-1, 0.), hadron_pt_sq(pTBin.size()-1, 0.); 
            std::vector <int> hadron_ct(pTBin.size()-1, 0), hadron_ct_sq(pTBin.size()-1, 0); 
            std::vector <double> hadron_phi_real(pTBin.size()-1, 0.), hadron_phi_imag(pTBin.size()-1, 0.), hadron_phi_real_sq(pTBin.size()-1, 0.), hadron_phi_imag_sq(pTBin.size()-1, 0.); 
        
            auto reader=make_shared<JetScapeReaderAscii>("/global/cscratch1/sd/td115/output/AuAu200/data_cut2_coef05/"+std::to_string(i_pT_hat_min)+"_i"+std::to_string(i_task)+".dat"); 
            std::ofstream hadron_output ("/global/cscratch1/sd/td115/output/AuAu200/observables/RAA/1e6_cut2_coef05/AuAu200_charged_hadron_"+std::to_string(i_pT_hat_min)+"_i"+std::to_string(i_task)+".txt");
            std::ofstream pt_output (("/global/cscratch1/sd/td115/output/AuAu200/observables/pTmean/1e6_cut2_coef05/AuAu200_charged_hadron_pT_"+std::to_string(i_pT_hat_min)+"_i"+std::to_string(i_task)+".txt").c_str()); 
            std::ofstream phi_output (("/global/cscratch1/sd/td115/output/AuAu200/observables/pt_phi/1e6_cut2_coef05/AuAu200_charged_hadron_pt_phi_"+std::to_string(i_pT_hat_min)+"_i"+std::to_string(i_task)+".txt").c_str()); 
/*
            double jetpTMin = 10., jetRadius = 0.2, jetAbsRapMax = 2., partonpTMin = 2.;
            double deltaRap = jetAbsRapMax * 2.; 
            fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, jetRadius); 
            vector <fjcore::PseudoJet> fjInputs;
            fjcore::Selector select_rapidity = fjcore::SelectorAbsRapMax(jetAbsRapMax);
            fjcore::Selector select_pt = fjcore::SelectorPtMin(partonpTMin); 
            fjcore::Selector select_both = select_rapidity && select_pt; 
*/
	    while (!reader->Finished())
	    {
                std::vector<int> hadron_ct_s(pTBin.size()-1, 0); 
                // std::vector<double> hadron_phi_real_s(pTBin.size()-1, 0.), hadron_phi_imag_s(pTBin.size()-1, 0.); 
                std::vector<double> hadron_pt_s(pTBin.size()-1, 0.); 
                reader->Next();
/*                fjInputs.resize(0);
                vector <fjcore::PseudoJet> inclusiveJets, sortedJets;
                fjcore::ClusterSequence clustSeq(reader->GetPartonsForFastJet(), jetDef);
                inclusiveJets = clustSeq.inclusive_jets(jetpTMin);
                vector <fjcore::PseudoJet> selected_jets = select_both(inclusiveJets);
                sortedJets = fjcore::sorted_by_pt(selected_jets);
                // std::cout << "parton size " << reader->GetPartonsForFastJet().size() << " jet size" << sortedJets.size() << "\n"; 
                for (unsigned int iJet = 0; iJet < sortedJets.size(); iJet++)
                {
                    double pT = sortedJets[iJet].perp(); 
                    if (pT < pTBin[0] || pT >= pTBin[pTBin.size()-1]) continue; 
                    double phi = sortedJets[iJet].phi(); 
                    int k = int((pT-10.)/2); 
		    hadron_ct_s[k]++; 
                    hadron_pt_s[k] += pT; 
                    phi_output << pT << " " << phi << "\n"; 
                }
               */ 
                auto hadrons = reader->GetHadrons(); 
	        auto pdghelper = JetScapeParticleBase::InternalHelperPythia.particleData;
                
                for (unsigned int iHadron = 0; iHadron < hadrons.size(); iHadron++)
                {
                    if (fabs(hadrons[iHadron]->eta()) > 1.) continue; 
                    auto ID = hadrons[iHadron]->pid(); 
                    auto charge = pdghelper.charge( ID ); 
                    if (charge == 0) continue; 
                    double pT = hadrons[iHadron]->pt();
                    if (pT < pTBin[0] || pT >= pTBin[pTBin.size()-1]) continue; 
                    auto phi = hadrons[iHadron]->phi();
                    int k = int((pT-10.)/2); 
                    // std::cout << pT << " " << ipT << "\n";  
                    // hadron_phi_real_s[k] += std::real(exp(2.*i*phi)); 
                    // hadron_phi_imag_s[k] += std::imag(exp(2.*i*phi)); 
		    hadron_ct_s[k]++; 
                    hadron_pt_s[k] += pT; 
                    phi_output << pT << " " << phi << "\n"; 
                }/*
                auto partons = reader->GetPartons(); 
                
                for (unsigned int iParton = 0; iParton < partons.size(); iParton++)
                {
                    double pT = partons[iParton]->pt();
                    if (pT < pTBin[0] || pT >= pTBin[pTBin.size()-1]) continue; 
                    auto phi = partons[iParton]->phi();
                    int k = int((pT-10.)/2); 
		    hadron_ct_s[k]++; 
                    hadron_pt_s[k] += pT; 
                    phi_output << pT << " " << phi << "\n"; 
                }*/
                for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
                {
                    hadron_ct[ipT] += hadron_ct_s[ipT]; 
                    // std::cout << hadron_ct_s[ipT] << "\n"; 
                    hadron_ct_sq[ipT] += hadron_ct_s[ipT] * hadron_ct_s[ipT]; 
                    // hadron_phi_real[ipT] += hadron_phi_real_s[ipT]; 
                    // hadron_phi_real_sq[ipT] += hadron_phi_real_s[ipT] * hadron_phi_real_s[ipT]; 
                    // hadron_phi_imag[ipT] += hadron_phi_imag_s[ipT]; 
                    // hadron_phi_imag_sq[ipT] += hadron_phi_imag_s[ipT] * hadron_phi_imag_s[ipT]; 
                    hadron_pt[ipT] += hadron_pt_s[ipT]; 
                    hadron_pt_sq[ipT] += hadron_pt_s[ipT] * hadron_pt_s[ipT]; 
                }
	    }
            for (unsigned int ipT = 0; ipT < pTBin.size() - 1; ipT++)
            {
                hadron_output << pTBin[ipT] << " " << hadron_ct[ipT] << " " << hadron_ct_sq[ipT] << "\n"; 
                // v2_output << pTBin[ipT] << " " << hadron_phi_real[ipT] << " " << hadron_phi_real_sq[ipT] << " " << hadron_phi_imag[ipT] << " " << hadron_phi_imag_sq[ipT] << "\n"; 
                pt_output << pTBin[ipT] << " " << hadron_pt[ipT] << " " << hadron_pt_sq[ipT] << "\n"; 
            }
	    reader->Close(); 
    cout<<"Finished!"<<endl; 
}
