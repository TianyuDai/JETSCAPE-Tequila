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
    double scale_list[1] = {1.}; 
    double alpha_list[1] = {0.3}; 
    double omegacut_list[1] = {1.}; 

    for (int ii = 0; ii < 1; ii++)
    {
        double scale = scale_list[ii]; 
        for (int jj = 0; jj < 1; jj++)
        {
            double alpha_s = alpha_list[jj]; 
            for (int kk = 0; kk < 1; kk++)
            {
                double omegacut = omegacut_list[kk]; 
                // auto reader=make_shared<JetScapeReaderAscii>("100GeV_gluon_muqperp"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"muomega"+std::to_string(omegacut)+"_QGP_inel_1fm.dat"); 
                auto reader=make_shared<JetScapeReaderAscii>("test_out.dat"); 
                // std::ofstream jet_output (("../../../../Result/JETSCAPE2.0/inel/general_test/qperp_test/100GeV_gluon_muqperp"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"muomega"+std::to_string(omegacut)+"_pureglue_elas_Lambda2_smallqperp.txt").c_str()); 
                // std::ofstream jet_output (("../../../../Result/JETSCAPE2.0/inel/cutoff_dependence_test/100GeV_gluon_muqperp"+std::to_string(scale)+"alpha"+std::to_string(alpha_s)+"muomega"+std::to_string(omegacut)+"_1fm_QGP_inel.txt").c_str()); 
                // std::ofstream jet_output ("../../../../Result/JETSCAPE2.0/general/Langevin/gluon_alpha_s0.3_p16_Langevin_200fm.txt"); 

                double total_energy = 0.; 
                int n = 0; 
                while (!reader->Finished())
                { 
                    reader->Next();
                    double energy = 0., px=0., py=0., pz=0.;
                    int id=21, stat=0;  
                    // jet_output<<"Event "<< reader->GetCurrentEvent()+1<<endl;
                    vector <shared_ptr <Parton> > partons;
                    partons = reader->GetPartons();
                    for (unsigned int i = 0; i < partons.size(); i++)
                   {
                       /*if (partons[i]->pstat() == 0) {*/energy = partons[i]->e();px = partons[i]->px(); py = partons[i]->py(); pz = partons[i]->pz(); id = partons[i]->pid(); 
                    stat = partons[i]->pstat();
                     total_energy += energy; 
                     n++; 
                    // jet_output << "energy: " << energy << ", px: " << px << ", py: " << py << ", id: " << id << ", stat: " << stat << endl;
// }
                        // if (abs(partons[i]->e()) > 1e-5) jet_output << partons[i]->e() << "\n";
                    // auto mShowers=reader->GetPartonShowers();
                    // for (int i=0;i<mShowers.size();i++)
                    //    for ( int ipart = 0; ipart< mShowers[i]->GetFinalPartons().size(); ++ipart)
                    //    {
                    //        Parton p = *mShowers[i]->GetFinalPartons().at(ipart);
                    //        if (p.e() > energy && p.pid() == 21) energy = p.e(); 
                    //      }
                    // if (partons[i]->pstat() != 1)
                    // if (partons[i]->pstat() == 1) std::cout << energy << "\n"; 
                           
                    } 
                }
                std::cout << "average energy is " << total_energy / 10000 << "\n"; 
                reader->Close(); 
            }
        }
    }
}
