#include "Tequila.h"
#include "IntTabulator.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include "tinyxml2.h"
#include "Srandom.h"
#include "FluidDynamics.h"
#include "TequilaMutex.h"
#include "JetScapeWriter.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp2d.h"
#include "gsl/gsl_spline2d.h"
#define nB(k) 1./(exp(k)-1.)
#define nF(k) 1./(exp(k)+1.)

#define MAGENTA "\033[35m"
#define hbarc 0.197327053

using namespace Jetscape;

Pythia8::Pythia Tequila::InternalHelperPythia ("IntentionallyEmpty",false);

// Register the module with the base class
RegisterJetScapeModule<Tequila> Tequila::reg("Tequila");


int Tequila::pLabelNew = 0; 
 
Tequila::Tequila()
{
    SetId("Tequila");
    VERBOSE(8);

    // create and set Tequila Mutex
    auto tequila_mutex = make_shared<TequilaMutex>(); 
    SetMutex(tequila_mutex); 
}

Tequila::~Tequila()
{
    VERBOSE(8);
}

void Tequila::Finish()
{
    free(za); 
    delete[] rate_ggg_p; 
    delete[] rate_gqq_p;
    delete[] rate_qqg_p;
    for(int ip=0;ip<nb_points_in_p; ip++)
    {
        for(int iomega=0;iomega<nb_points_in_omega; iomega++) 
        {
            delete[] differential_rate_ggg_p_omega_qperp[ip][iomega]; 
            delete[] differential_rate_gqq_p_omega_qperp[ip][iomega];
            delete[] differential_rate_qqg_p_omega_qperp[ip][iomega];
        }
        delete[] differential_rate_ggg_p_omega_qperp[ip];
        delete[] differential_rate_gqq_p_omega_qperp[ip];
        delete[] differential_rate_qqg_p_omega_qperp[ip];
    }
    delete[] differential_rate_ggg_p_omega_qperp;
    delete[] differential_rate_gqq_p_omega_qperp;
    delete[] differential_rate_qqg_p_omega_qperp;
}

void Tequila::Init()
{
    JSINFO<<"Intialize Tequila ...";

    double deltaT = 0.; 
    deltaT = GetXMLElementDouble({"Eloss", "deltaT"}); 

    string s = GetXMLElementText({"Eloss", "Tequila", "name"});
    JSDEBUG << s << " to be initilizied ...";

    Q0 = GetXMLElementDouble({"Eloss", "Tequila", "Q0"});
    alphas_soft = GetXMLElementDouble({"Eloss", "Tequila", "alphas_soft"});
    alphas_hard_elas = GetXMLElementDouble({"Eloss", "Tequila", "alphas_hard_elas"});
    alphas_hard_inel = GetXMLElementDouble({"Eloss", "Tequila", "alphas_hard_inel"});
    pcut = GetXMLElementDouble({"Eloss", "Tequila", "pcut"});
    hydro_Tc = GetXMLElementDouble({"Eloss", "Tequila", "hydro_Tc"});
    recoil_on = GetXMLElementInt({"Eloss", "Tequila", "recoil_on"});

    muqperp_over_T = GetXMLElementDouble({"Eloss", "Tequila", "muqperp_over_T"});
    muomega_over_T = GetXMLElementDouble({"Eloss", "Tequila", "muomega_over_T"});
    // qhat_coef = GetXMLElementDouble({"Eloss", "Tequila", "qhat_coef"});

    // Path to additional data
    path_to_tables = GetXMLElementText({"Eloss", "Tequila", "path"});

    g_soft = sqrt(4.*M_PI*alphas_soft); 
    g_hard_elas = sqrt(4.*M_PI*alphas_hard_elas); 
    g_hard_inel = sqrt(4.*M_PI*alphas_hard_inel); 

    const double alpha_EM=1./137.; //Not currently a parameter
    hydro_tStart = 0.4; 
    ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 }; 
	
    // Load elastic rate
    JSINFO << "Load elastic rate... "; 
    LoadElasticTables(); 

    // Load differential rate
    allocate_memory_for_radiative_rate_table();

    JSINFO << "Loading differential collinear rate...\n";
    // const std::string location_of_pretabulated_collinear_rates="./Tequila/";
    // Either compute or read from file the collinear differential rates and load them into member array "differential_rate_p_omega_qperp[][][]"
    load_differential_rate(alphas_hard_inel, alpha_EM, nf, path_to_tables);
    // Compute total rate from differential rate and save in member array "rate_p[]"
    JSINFO << "Computing integrated collinear rate...\n";
    evaluate_integrated_rate(muomega_over_T, differential_rate_qqg_p_omega_qperp,rate_qqg_p, qqg); 
    evaluate_integrated_rate(muomega_over_T, differential_rate_ggg_p_omega_qperp,rate_ggg_p, ggg);
    evaluate_integrated_rate(muomega_over_T, differential_rate_gqq_p_omega_qperp,rate_gqq_p, gqq);
    JSINFO << "Done computing integrated collinear rate.\n";

}

void Tequila::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

    VERBOSESHOWER(5)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<< pIn.size();
    int Id, IdNew=21;
    int pStat, pStatNew, pStatTmp=-1; 		   // status of parton: 
			   // 0: normal parton, 1: recoil parton or radiated parton, 
			   // -1: sampled thermal parton (negative)
    int pLabel;     // particle number
    double pAbs, px, py, pz;   // momentum of initial parton (pIn)
    double pRest, pxRest, pyRest, pzRest;      // momentum in the rest frame of fluid cell (pIn)
    double pRestNew; 
    FourVector pVec, pVecRest; 
    FourVector pVecRestNew, pVecRestNewest, pVecNewest; 	// 4 vector momentum of parton after hard interaction & diffusion
    double k, kRest=0.;           // momentum of radiated parton (pOut)
    FourVector kVec, kVecRest, kVecRestNew, kVecRestNewest, kVecNewest;
    double pThermal, pxThermal;   // momentum of thermal parton (pOut)
    double pyThermal, pzThermal;
    FourVector pVecThermal, pVecThermalRest; 
    double pRecoil, pxRecoil;     // momentum of recoil parton (pOut)
    double pyRecoil, pzRecoil;    
    double pRecoilRest = 0.;
    FourVector pVecRecoilNewest, pVecRecoilRest, pVecRecoilRestNewest;
    FourVector pVecLangevinTmp, kVecLangevinTmp, pVecRecoilLangevinTmp;;  

    double xx, yy, zz, tt;         // position of initial parton (pIn)
    FourVector xVec;           // 4 vector for position (for next time step!)

    double eta;                // pseudo-rapidity
    double velocity_jet[4];    // jet velocity for MATTER
	
    // flow info
    double vx, vy, vz;         // 3 components of flow velocity
    double T;                  // Temperature of fluid cell
    double beta, gamma; 

    double omega = 0., qperp = 0., q = 0.; // transfer energy and transverse momentum, 
	                           // qperp is actually \tilde{q_\perp}=sqrt(q^2 - omega^2) here, we should convert it to the real q_\perp
    FourVector qVec; 
    for (int i=0;i<pIn.size();i++)
    {
        // Reject photons
        if (pIn[i].pid()==photonid)
            continue;
        if (abs(pIn[i].pid()) > 3 && pIn[i].pid() != 21)
            continue;
 
        Id = pIn[i].pid();
        pStat = pIn[i].pstat(); 
        // do nothing for negative particles
        if (pStat < 0) continue;
 
        px = pIn[i].px();
        py = pIn[i].py();
        pz = pIn[i].pz();

        // In Tequila, particles are all massless and on-shell
        M = InternalHelperPythia.particleData.m0( Id ); 
        pAbs = sqrt(px*px+py*py+pz*pz+M*M);
        pVec = FourVector ( px, py, pz, pAbs );

        tt = pIn[i].x_in().t(); 
        xx = pIn[i].x_in().x() + (Time-tt)*px/pAbs;
        yy = pIn[i].x_in().y() + (Time-tt)*py/pAbs;
        zz = pIn[i].x_in().z() + (Time-tt)*pz/pAbs;
        eta = pIn[i].eta();

        // Extract fluid properties
        std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

        GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
        VERBOSE(8) << MAGENTA << "Temperature from Brick (Signal) = " << check_fluid_info_ptr->temperature;

        vx = check_fluid_info_ptr->vx;
        vy = check_fluid_info_ptr->vy;
        vz = check_fluid_info_ptr->vz;
        T = check_fluid_info_ptr->temperature;
        static Pythia8::Pythia InternalHelperPythia;
        // Only accept low t particles
        if (pIn[i].t() > Q0*Q0 + rounding_error || T < hydro_Tc) continue; 
        TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.
        
        pLabel = pIn[i].plabel();
        if(pLabel == 0) 
        {
            IncrementpLabel();
            pIn[i].set_label(pLabelNew);
            pLabel = pLabelNew;
        }

        beta = sqrt( vx*vx + vy*vy + vz*vz );
        gamma = 1./sqrt(1.-beta*beta); 

        // Set momentum in fluid cell's frame
        fourvec pIn_4vec = fourvec{pAbs, px, py, pz}; 
        fourvec pIn_cell = pIn_4vec.boost_to(vx, vy, vz); 
        pVecRest = FourVector(pIn_cell.x(), pIn_cell.y(), pIn_cell.z(), pIn_cell.t()); 
        pRest = pVecRest.t(); 
        pxRest = pVecRest.x(); 
        pyRest = pVecRest.y(); 
        pzRest = pVecRest.z(); 

        if (pRest < eLossCut) continue; 

        VERBOSE(8) << MAGENTA
            << "Time = " << Time << " Id = " << Id << " T = " << T
            << " pAbs = " << pAbs << " " << px << " " << py << " " << pz 
            << " | position = " << xx << " " << yy << " " << zz;
      	
      	// Pass the address of deltaT to DetermineProcess, if deltaT is changed, the changed value could be passed back. Calculate xVec after passing back the real deltaT. 
      	// time step should be in the rest frame. 
        double deltaTRest = deltaT / gamma; 
        Lambda = std::min({pRest, 2.*T*elas_omega_over_T_pos_max, 2.*sqrt(3.*pRest*T), 2.*pcut});

        process_type process = DetermineProcess(pRest, T, deltaTRest, Id);

      	xVec = FourVector( xx+px/pAbs*deltaT, yy+py/pAbs*deltaT, zz+pz/pAbs*deltaT, Time+deltaT );
        velocity_jet[0]=1.0;
        velocity_jet[1]=pIn[i].jet_v().x();
        velocity_jet[2]=pIn[i].jet_v().y();
        velocity_jet[3]=pIn[i].jet_v().z();
        
        IntTabulator inttabulator;
        if (process == none)
      	    pVecRestNew = pVecRest;
      	else if (process == gg || process == gq || process == qg || process == qq || process == qqp || process == qqb)
      	{
      	    omega = Energy_Transfer(pRest, T, process); 
	    qperp = TransverseMomentum_Transfer(pRest, omega, T, process);
	    // qperp = 0.1;
            pVecRestNew = Momentum_Update(omega, qperp, T, pVecRest); 
            pRestNew = pVecRestNew.t();
            /* 
            if (recoil_on)
            {
                qVec = FourVector(pVecRest.x()-pVecRestNew.x(), pVecRest.y()-pVecRestNew.y(), pVecRest.z()-pVecRestNew.z(), pVecRest.t()-pVecRestNew.t()); 
                pVecThermalRest = getThermalVec(qVec, T, Id); 
                pVecRecoilRest = FourVector(qVec.x()+pVecThermalRest.x(), qVec.y()+pVecThermalRest.y(), qVec.z()+pVecThermalRest.z(), qVec.t()+pVecThermalRest.t()); 
                pRecoilRest = pVecRecoilRest.t();
                if (process == gg || process == qg) IdNew = 21; 
                else if (process == qq) IdNew = Id; 
                else if (process == gq)
                {
                    double r = ZeroOneDistribution(*GetMt19937Generator());
                    if (r < 1./6.) IdNew = 1;
                    else if (r < 2./6.) IdNew = 2;
                    else if (r < 3./6.) IdNew = 3;
                    else if (r < 4./6.) IdNew = -1;
                    else if (r < 5./6.) IdNew = -2;
                    else IdNew = -3;    
                }
                else if (process == qqb) IdNew = -1 * Id; 
                else if (process == qqp)
                {
                    double r = ZeroOneDistribution(*GetMt19937Generator()); 
                    if (r < 1./4.) IdNew = (abs(Id) + 1) % 3 + 1; 
                    else if (r < 2./4.) IdNew = (abs(Id) + 2) % 3 + 1; 
                    else if (r < 3./4.) IdNew = -1 * ((abs(Id) + 1) % 3 + 1); 
                    else IdNew = -1 * ((abs(Id) + 2) % 3 + 1); 
                }
                pStatNew = 1; 
            }*/
        }
      	else if (process == GqQg || process == GgQbq)
      	{
      	    omega = Energy_Transfer(pRest, T, process)*T; 
	        qperp = 0.; 
	        pVecRestNew = Momentum_Update(omega, qperp, T, pVecRest); 
	        // choose the Id of new qqbar pair. Note that we only deal with nf = 3
	        double r = ZeroOneDistribution(*GetMt19937Generator());
	        if (r < 1./6.) Id = 1;
	        else if (r < 2./6.) Id = 2;
	        else if (r < 3./6.) Id = 3;
	        else if (r < 4./6.) Id = -1; 
	        else if (r < 5./6.) Id = -2; 
	        else Id = -3; 
	    }
      	else if (process == QgGq || process == QbqGg)
      	{
      	    omega = Energy_Transfer(pRest, T, process)*T; 
	    qperp = 0.; 
	    pVecRestNew = Momentum_Update(omega, qperp, T, pVecRest);
            Id = 21;  
	}
        else if (process == gg_split || process == gq_split || process == qg_split || process == qq_split || process == qqp_split || process == qqb_split)
        {
            double x = xSampling(pRest, T, process); 
            omega = pRest * x; 
            // qperp = 0.1; 
	        qperp = TransverseMomentum_Transfer_Split(pRest, omega, T, process);
            pVecRestNew = Momentum_Update(omega, qperp, T, pVecRest);
            if (recoil_on)
            {
                qVec = FourVector(pVecRest.x()-pVecRestNew.x(), pVecRest.y()-pVecRestNew.y(), pVecRest.z()-pVecRestNew.z(), pVecRest.t()-pVecRestNew.t());
                pVecThermalRest = getThermalVec(qVec, T, Id);
                pVecRecoilRest = FourVector(qVec.x()+pVecThermalRest.x(), qVec.y()+pVecThermalRest.y(), qVec.z()+pVecThermalRest.z(), qVec.t()+pVecThermalRest.t()); 
                pRecoilRest = pVecRecoilRest.t(); 
             
                if (process == gg_split || process == qg_split) IdNew = 21; 
                else if (process == gq_split)
                {
                    double r = ZeroOneDistribution(*GetMt19937Generator());
                    if (r < 1./6.) IdNew = 1;
                    else if (r < 2./6.) IdNew = 2;
                    else if (r < 3./6.) IdNew = 3;
                    else if (r < 4./6.) IdNew = -1;
                    else if (r < 3./6.) IdNew = 3;
                    else if (r < 4./6.) IdNew = -1;
                    else if (r < 5./6.) IdNew = -2;
                    else IdNew = -3;
                }
                else if (process == qq_split)
                    IdNew = Id; 
                else if (process == qqp_split)
                {
                    double r = ZeroOneDistribution(*GetMt19937Generator());
                    if (r < 1./4.) IdNew = (abs(Id)+1) % 3 + 1;
                    else if (r < 2./4.) IdNew = (abs(Id)+2) % 3 + 1;
                    else if (r < 3./4.) IdNew = -1*((abs(Id)+1) % 3 + 1);
                    else IdNew = -1*((abs(Id)+2) % 3 + 1);
                }
                else if (process == qqb_split)
                    IdNew = -1*Id; 
                pStatNew = 1; 
            }
        }
        else if (process == ggg)
      	{ 
            if (pRest/T < AMY_p_over_T_cut) return;
            // sample radiated parton's momentum
            sample_dgamma_dwdq(pRest, T,differential_rate_ggg_p_omega_qperp, omega, qperp, ggg);
            kRest = omega;
            if(kRest > pRest) return;

            // final state parton's momentum
            pRestNew = pRest - kRest;

            pVecRestNew.Set( (pxRest/pRest)*pRestNew, (pyRest/pRest)*pRestNew, (pzRest/pRest)*pRestNew, pRestNew );

            kVecRest.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest );
            IdNew = 21;
            pStatNew = 1;  
        }
        else if (process == gqq)
        {
            if (pRest/T < AMY_p_over_T_cut) return;

            // sample radiated parton's momentum
            sample_dgamma_dwdq(pRest, T,differential_rate_gqq_p_omega_qperp, omega, qperp, gqq); 
            kRest = omega; 
            if(kRest > pRest) return;

            // final state parton's momentum
            pRestNew = pRest - kRest;

            // choose the Id of new qqbar pair. Note that we only deal with nf = 3
            double r = ZeroOneDistribution(*GetMt19937Generator());
            if (r < 1./6.) Id = 1;
            else if (r < 2./6.) Id = 2;
            else if (r < 3./6.) Id = 3;
            else if (r < 4./6.) Id = -1; 
            else if (r < 5./6.) Id = -2; 
            else Id = -3; 

            pVecRestNew.Set( (pxRest/pRest)*pRestNew, (pyRest/pRest)*pRestNew, (pzRest/pRest)*pRestNew, pRestNew );

            kVecRest.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest );
            IdNew = -1 * Id;
            pStatNew = 1;  
        }
        else if (process == qqg)
        {
            if (pRest/T < AMY_p_over_T_cut) return;

            // sample radiated parton's momentum
            sample_dgamma_dwdq(pRest, T,differential_rate_qqg_p_omega_qperp, omega, qperp, qqg); 
            kRest = omega; 
            if(kRest > pRest) return;

            // final state parton's momentum
            pRestNew = pRest - kRest;

            pVecRestNew.Set( (pxRest/pRest)*pRestNew, (pyRest/pRest)*pRestNew, (pzRest/pRest)*pRestNew, pRestNew );

            kVecRest.Set( (pxRest/pRest)*kRest, (pyRest/pRest)*kRest, (pzRest/pRest)*kRest, kRest );
            IdNew = 21; 
            pStatNew = 1; 
        }
        else pVecRestNew = pVecRest;
        // pVecRestNew = pVecRest; 
        
        // diffusion process
        pVecRestNewest = Langevin_Update(deltaTRest / hbarc, T, pVecRestNew, Id);

        // pVecRestNewest = pVecRestNew;
        if (pVecRestNewest.t() > pcut)
        {
            fourvec pOut_4vec, kOut_4vec; 
            // boost the updated parton from the rest frame back to the fluid frame
            pOut_4vec = fourvec{pVecRestNewest.t(), pVecRestNewest.x(), pVecRestNewest.y(), pVecRestNewest.z()}; 
            pOut_4vec = pOut_4vec.boost_back(vx, vy, vz); 
            pVecNewest = FourVector(pOut_4vec.x(), pOut_4vec.y(), pOut_4vec.z(), pOut_4vec.t());
 
            // push the updated parton back into the parton list
            pOut.push_back(Parton(pLabel, Id, pStat, pVecNewest, xVec)); 
            pOut[pOut.size()-1].set_form_time(0.);
            pOut[pOut.size()-1].set_jet_v(velocity_jet); 
        }
        if (kRest > epsilon)
        {
            kVecRestNewest = Langevin_Update(deltaTRest / hbarc, T, kVecRest, IdNew); 
            // kVecRestNewest = kVecRest; 

            fourvec kOut_4vec; 
            kOut_4vec = fourvec{kVecRestNewest.t(), kVecRestNewest.x(), kVecRestNewest.y(), kVecRestNewest.z()}; 
            kOut_4vec = kOut_4vec.boost_back(vx, vy, vz);

            kVecNewest = FourVector(kOut_4vec.x(), kOut_4vec.y(), kOut_4vec.z(), kOut_4vec.t());
            if (kVecNewest.t() > pcut)
            { 
                IncrementpLabel();
                pOut.push_back(Parton(pLabelNew, IdNew, 1, kVecNewest, xVec));
                pOut[pOut.size()-1].set_form_time(0.);
                pOut[pOut.size()-1].set_jet_v(velocity_jet);
            }
       }
       
        // push the recoil parton back into the parton list
       if (pRecoilRest > epsilon) 
       {
            pVecRecoilRestNewest = Langevin_Update(deltaTRest / hbarc, T, pVecRecoilRest, IdNew); 
            // pVecRecoilRestNewest = pVecRecoilRest; 

            fourvec pRecoilOut_4vec; 
            pRecoilOut_4vec = fourvec{pVecRecoilRestNewest.t(), pVecRecoilRestNewest.x(), pVecRecoilRestNewest.y(), pVecRecoilRestNewest.z()}; 
            pRecoilOut_4vec = pRecoilOut_4vec.boost_back(vx, vy, vz);

            pVecRecoilNewest = FourVector(pRecoilOut_4vec.x(), pRecoilOut_4vec.y(), pRecoilOut_4vec.z(), pRecoilOut_4vec.t()); 
            // std::cout << "recoil rest " << pVecRecoilRest.x() << " " << pVecRecoilRest.x() << " " << pVecRecoilRest.z() << " " << pVecRecoilRest.t() << "\n"; 
            // std::cout << "recoil lab " << pVecRecoilNewest.x() << " " << pVecRecoilNewest.y() << " " << pVecRecoilNewest.z() << " " << pVecRecoilNewest.t() << "\n"; 
            if (pVecRecoilNewest.t() > pcut)
            {
                // std::cout << pVecRecoilNewest.t() << "\n"; 
                IncrementpLabel(); 
                pOut.push_back(Parton(pLabelNew, IdNew, 1, pVecRecoilNewest, xVec)); 
                pOut[pOut.size()-1].set_form_time(0.); 
                pOut[pOut.size()-1].set_jet_v(velocity_jet);
            }
 
        }
      	return; 
    }
}

process_type Tequila::DetermineProcess(double pRest, double T, double deltaTRest, int Id)
{
    double dt = deltaTRest / hbarc; 
    double rateTotal; 
    double rate[nTotalProcess] = {0.};
 
    // Elastic rate	
    for (int i = gg; i <= qqb; i++)
    {
        process_type process = static_cast<process_type>(i); 
        rate[i] = interpolatorElasticTotalRate(elas_omega_over_T_pos_max*T, T, process); 
        rate[i] += casimir[i] * pow(g_hard_elas*g_hard_elas*T, 2) / (2*elas_omega_over_T_pos_max*T);
        rate[i] -= casimir[i] * pow(g_hard_elas*g_hard_elas*T, 2) / Lambda;     // Add the cutoff Lambda
        rate[i] -= interpolatorElasticEnergyLoss(elas_omega_over_T_pos_max*T, T, process) / pRest;
        rate[i] -= 1. / 96 / M_PI / pRest * pow(g_hard_elas*g_hard_elas*T, 2) * Toverp_c_ln[i]*log(elas_omega_over_T_pos_max);
        rate[i] += 1. / 96 / M_PI / pRest * pow(g_hard_elas*g_hard_elas*T, 2) * Toverp_c_ln[i]*log(Lambda/T/2);
        // rate[i] += casimir[i] * pow(g*g*T, 2) / (2*elas_omega_over_T_pos_max*T); 
        // rate[i] -= casimir[i] * pow(g*g*T, 2) / Lambda; 
    }
    // JSINFO << "rate gg " << rate[gg] << " rate gq " << rate[gq] << " rate qg " << rate[qg] << " rate qq " << rate[qq] << " rate qqb " << rate[qqb]; 
    // std::cout << "gg " << rate[gg] << " gq " << rate[gq] << " qg " << rate[qg] << " qqp " << rate[qqp] << " qqb " << rate[qqb] << " GqQg " << rate[GqQg] << " QgGq " << rate[QgGq] << " GgQbq " << rate[GgQbq] << " QbqGg " << rate[QbqGg] << "\n"; 

    for (int i = GqQg; i <= QbqGg; i++)
    {
        process_type process = static_cast<process_type>(i); 
        rate[i] = interpolatorElasticTotalRate(elas_omega_over_T_pos_max*T, T, process) * T / pRest;
        // JSINFO << "conversion " << rate[i] << "\n"; 
        rate[i] += 1. / 96 / M_PI * pow(g_hard_elas, 4) * T * Toverp_c_ln[i]*log(elas_omega_over_T_pos_max) * T /pRest;
        rate[i] -= 1. / 96 / M_PI * pow(g_hard_elas, 4) * T * Toverp_c_ln[i]*log(Lambda/2/T) * T / pRest;
        // JSINFO << "coef " << g << " " << Toverp_c_ln[i] << " " << Lambda << " " << log(elas_omega_over_T_pos_max) << "\n"; 
    }
    // JSINFO << "rate GqQg " << rate[GqQg] << " QbqGg " << rate[QbqGg] << "\n"; 

    for (int i = gg_split; i <= qqb_split; i++)
    {	
        rate[i] = splitting_c_lambda[i - gg_split]*pRest/Lambda + splitting_c_p[i - gg_split] + splitting_c_ln[i - gg_split]*log(2*pRest/Lambda); 
        if (i == gg_split || i == qg_split)
            rate[i] *= pow(g_hard_elas, 4)*pow(T, 2)/96/M_PI/pRest;
        else
            rate[i] *= pow(g_hard_elas, 4)*pow(T, 2)/1536/M_PI/pRest;
        // rate[i] *= pow(g, 4)*pow(T, 2)/M_PI/pRest; 
    }

    // JSINFO << "gg_split " << rate[gg_split] << " " << "gq_split " << rate[gq_split] << " " << "qg_split " << rate[qg_split] << " " << "qqp_split " << rate[qqp_split] << " " << "qqb_split " << rate[qqb_split] << " qq_split " << rate[qq_split] << "\n"; 
/*    JSINFO << BOLDGREEN << "pRest " << pRest << " Lambda " << Lambda << " gg_split " << rate[gg_split] << " gg_elas " << rate[gg] << "\n"; 
    JSINFO << BOLDGREEN << "pRest " << pRest << " Lambda " << Lambda << " gq_split " << rate[gq_split] << " gq_elas " << rate[gq] << "\n"; 
    JSINFO << BOLDGREEN << "pRest " << pRest << " Lambda " << Lambda << " qqp_split " << rate[qqp_split] << " qqp_elas " << rate[qqp] << "\n"; 
    JSINFO << BOLDGREEN << "pRest " << pRest << " Lambda " << Lambda << " qg_split " << rate[qg_split] << " qg_elas " << rate[qg] << "\n"; 
*/
    // Inelastic rate
    rate[qqg] = rate_inel(pRest, T, rate_qqg_p); 
    rate[gqq] = rate_inel(pRest, T, rate_gqq_p); 
    rate[ggg] = rate_inel(pRest, T, rate_ggg_p); 

    // JSINFO << "rate qqg " << rate[qqg] << " rate gqq " << rate[gqq] << " rate ggg " << rate[ggg]; 

    if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
    {
        double totalQuarkProb = 0.; 
        if (pRest/T > AMY_p_over_T_cut)
            totalQuarkProb += rate[qqg]*dt;
        totalQuarkProb += (rate[qg] + rate[qq] + rate[qqp] + rate[qqb] + rate[qg_split] + rate[qq_split] + rate[qqp_split] + rate[qqb_split]) * dt; 
        // warn if total probability exceeds 0.2
        if (totalQuarkProb > 1.)
            JSWARN << " : Total Probability for quark processes exceeds 1 (" << totalQuarkProb << "). " << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01)."; 	
	
        double accumProb = 0.; 
        double nextProb = 0.; 
        double Prob = 0.; 
        double randProb = ZeroOneDistribution(*GetMt19937Generator()); 

        if (randProb < totalQuarkProb)
        {
            Prob = rate[qg] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qg; 

            accumProb += Prob; 
            Prob = rate[qq] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qq;

            accumProb += Prob; 
            Prob = rate[qqp] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqp; 

            accumProb += Prob; 
            Prob = rate[qqb] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqb;

            accumProb += Prob;
            Prob = rate[QgGq] * dt;
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return QgGq;

            accumProb += Prob;
            Prob = rate[QbqGg] * dt;
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return QbqGg;

            accumProb += Prob; 
            Prob = rate[qg_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qg_split; 

            accumProb += Prob; 
            Prob = rate[qq_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qq_split; 

            accumProb += Prob; 
            Prob = rate[qqp_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqp_split; 

            accumProb += Prob; 
            Prob = rate[qqb_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqb_split;
            /* 
            accumProb += Prob; 
            Prob = rate[qqbp_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return qqbp_split; 
            */
            accumProb += Prob; 
            Prob = rate[qqg] * dt; 
            if (pRest/T > AMY_p_over_T_cut && accumProb <= randProb && randProb < (accumProb + Prob)) return qqg; 
        }
        else
            return none; 
    }
    else if (Id == 21)
    {
        double totalGluonProb = 0.; 
        if (pRest/T > AMY_p_over_T_cut) 
            totalGluonProb += (rate[gqq] + rate[ggg])*dt;
        totalGluonProb += (rate[gg] + rate[gq] + rate[gg_split] + rate[gq_split]) * dt; 
        if (totalGluonProb > 1.)
            JSWARN << " : Total Probability for gluon processes exceeds 1 (" << totalGluonProb << "). " << " : Most likely this means you should choose a smaller deltaT in the xml (e.g. 0.01)."; 

        double accumProb = 0.; 
        double nextProb = 0.; 
        double Prob = 0.; 
        double randProb = ZeroOneDistribution(*GetMt19937Generator()); 

        if (randProb < totalGluonProb)
        {
            Prob = rate[gg] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return gg; 

            accumProb += Prob; 
            Prob = rate[gq] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return gq;

            accumProb += Prob;
            Prob = rate[GqQg] * dt;
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return GqQg;

            accumProb += Prob;
            Prob = rate[GgQbq] * dt;
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return GgQbq;

            accumProb += Prob; 
            Prob = rate[gg_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return gg_split; 

            accumProb += Prob; 
            Prob = rate[gq_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return gq_split;
            /*
            accumProb += Prob; 
            Prob = rate[ggqqb_split] * dt; 
            if (accumProb <= randProb && randProb < (accumProb + Prob)) return ggqqb_split;
            */
            accumProb += Prob; 
            Prob = rate[ggg] * dt; 
            if (pRest/T > AMY_p_over_T_cut && accumProb <= randProb && randProb < (accumProb + Prob)) return ggg; 

            accumProb += Prob; 
            Prob = rate[gqq] * dt; 
            if (pRest/T > AMY_p_over_T_cut && accumProb <= randProb && randProb < (accumProb + Prob)) return gqq; 

        }
        else return none; 
    }
    return none; 
}

FourVector Tequila::Momentum_Update(double omega, double qperp, double T, FourVector pVec)
{
    double p = pVec.t(); 
    double pp = p - omega; 
    double sintheta = fabs(qperp/pp); 
    double phi = 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator());
    fourvec pOut, pIn; 
    FourVector pVecNew; 
    pIn = fourvec{pVec.t(), pVec.x(), pVec.y(), pVec.z()}; 
    // Assume pVec is along z axis. 
    pOut = fourvec{fabs(pp), fabs(pp)*sintheta*cos(phi), fabs(pp)*sintheta*sin(phi), pp*sqrt(1-sintheta*sintheta)}; 
    pOut = pOut.rotate_back(pIn); 
    pVecNew = FourVector(pOut.x(), pOut.y(), pOut.z(), pOut.t()); 
    return pVecNew; 
}


FourVector Tequila::getThermalVec(FourVector qVec, double T, int id)
{
    int kind = 1;               // 1: fermion, -1: boson
    if (fabs(id) < 4)  kind = 1;
    else if (id == 21) kind = -1;

    FourVector kVec;            // new thermal parton's momentum
    double k, k_min;            // (minimum) mangitude of thermal momentum
    double cosTheta, sinTheta;  // angle between vector q and k

    double q = sqrt(qVec.x()*qVec.x()+qVec.y()*qVec.y()+qVec.z()*qVec.z());
    double omega = qVec.t();

    // if momentum transfer is on-shell or time-like
    // return null FourVector
    if (q - fabs(omega) < 1e-5)
        return FourVector( 0., 0., 0., 0.);

    // minimum momenum of thermal parton k that makes recoil parton on-shell
    k_min = (q-omega)/2.;

    // sampled magnitude of thermal momentum
    // JSINFO << "kmin " << k_min << " kind " << kind; 
    k = getThermal(k_min, T, kind);
    // JSINFO << "thermal energy is " << k; 
    // Theta is the angle between k and q
    cosTheta = -1. * (2.*k*omega - q*q + omega*omega)/(2.*k*q);
    sinTheta = sqrt(1.-cosTheta*cosTheta);
    double phi = 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator()); 

    fourvec qVec_fourvec, kVec_fourvec; 
    qVec_fourvec = fourvec{qVec.t(), qVec.x(), qVec.y(), qVec.z()}; 
    // Assume qVec is along z axis. 
    kVec_fourvec = fourvec{fabs(k), fabs(k)*sinTheta*cos(phi), fabs(k)*sinTheta*sin(phi), k*sqrt(1-sinTheta*sinTheta)}; 
    kVec_fourvec = kVec_fourvec.rotate_back(qVec_fourvec); 

    kVec = FourVector(kVec_fourvec.x(), kVec_fourvec.y(), kVec_fourvec.z(), kVec_fourvec.t()); 
    return kVec; 
}

double Tequila::getThermal(double k_min, double T, int kind)
{
    // this algorithm uses 5.5/(1+px2) as envelope function and then uses the rejection method
    // tan (Pi/2*x) is the inverse of the integral of the envelope function
    // this can be improved
    // kind: -1 = Boson
    //        1 = Fermion
    double p;
    double gm = 5.5;
    bool repeat = true;
    double gx = 0.;
    double px;
    double px2;
    double ex;
    int n_try = 0, max_try = 10000;
    if (k_min > 3.) return k_min;  
    do 
    {
        px = k_min/T + tan (M_PI * ZeroOneDistribution(*GetMt19937Generator()) / 2.0);
        px2 = px*px;
        ex = sqrt (px2);
        if (kind == 0) gx = px2 * (1.0 + px2) * exp (- ex);
        else if (kind == -1) gx = px2 * (1.0 + px2) * 1/(exp(ex)-1);
        else if (kind == +1) gx = px2 * (1.0 + px2) * 1/(exp(ex)+1);
        if ( ZeroOneDistribution(*GetMt19937Generator()) < gx / gm)
        {
            repeat = false;
            p = ex;
        } 
        n_try++; 
        if (n_try > max_try)
        {
            JSWARN << "Cannot find a proper thermal energy! kmin is " << k_min; 
            return k_min; 
        }
    } while (repeat); 

    return p*T;  // p*T is in [GeV]
}

double Tequila::qhatpara(double E, double T, int id)
{
    double CR; 
    if (id == 21) CR = CA; 
    else if (std::abs(id) == 1 || std::abs(id) == 2 || std::abs(id) == 3) CR = CF; 
    else {JSWARN << "Strange particle! ID = " << id; CR = CF; }
    double mD = sqrt(std::pow(g_soft*T, 2)*(nc/3. + nf/6.)); 
    double Minf = sqrt(pow(mD, 2)/2.); 
    double qpara_elas = std::pow(g_soft*Minf, 2)*CR*T/(2.*M_PI)*log(1.+pow(muqperp_over_T*T/Minf, 2))/2.; 
    double qpara_inel = std::pow(g_soft, 4)*CR*CA*std::pow(T, 3)*muomega_over_T*(2-ln2)/(4.*std::pow(M_PI, 3)); 
    return qpara_elas + qpara_inel; 
}

double Tequila::qhatperp(double E, double T, int id)
{
    double CR; 
    if (id == 21) CR = CA; 
    else if (std::abs(id) == 1 || std::abs(id) == 2 || std::abs(id) == 3) CR = CF; 
    else {JSWARN << "Strange particle! ID = " << id; CR = CF; }
    double mD = sqrt(std::pow(g_soft*T, 2)*(nc/3. + nf/6.));
    double qperp_elas = std::pow(g_soft*mD, 2) * CR * T / (2.*M_PI) * log(1.+pow(muqperp_over_T*T/mD, 2))/2.; 
    double qperp_inel = 0.;  
    return qperp_elas + qperp_inel;
}

FourVector Tequila::Langevin_Update(double dt, double T, FourVector pIn, int id)
{
    // imaging rotating to a frame where pIn lies on z-axis
    double E0 = pIn.t();
    double p0 = std::sqrt(E0*E0 - M*M + 1e-9);
    double kt = qhatperp(E0, T, id) / 2.;
    double kl = qhatpara(E0, T, id); 
    double drag = kl/(2.*p0*T)+1./(2.*p0*p0)*(kt*2.-kl*2.); //eta_D
		   
    double white_noise_holder[3];
    for (size_t i=0; i<3; ++i) 
        white_noise_holder[i] = Srandom::white_noise(Srandom::gen);

    double Ct = std::sqrt(kt*dt);
    double Cl = std::sqrt(kl*dt);
	
    fourvec pOut; 
    pOut.a[1] = Ct * white_noise_holder[0];
    pOut.a[2] = Ct * white_noise_holder[1];
    pOut.a[3] = p0 * (1. - drag * dt) + Cl * white_noise_holder[2];
    pOut.a[0] = std::sqrt(M*M + std::pow(pOut.x(),2) + std::pow(pOut.y(),2) + std::pow(pOut.z(),2) );

    // rotate back to the original frame
    pOut = pOut.rotate_back(fourvec{pIn.t(), pIn.x(), pIn.y(), pIn.z()});
    return FourVector(pOut.x(), pOut.y(), pOut.z(), pOut.t()); 
}

double Tequila::TransverseMomentum_Transfer(double pRest, double omega, double T, process_type process)
{
    double pRest_over_T = pRest / T; 
    double omega_over_T = omega / T;
    // JSINFO << pRest_over_T << " " << omega_over_T;  
    // change the sampling range to make q_perp^2 dGamma/domega/dqperp2 to be a fair value
    double limitMax = std::min(pow(std::min(elas_qperp_over_T_max, fabs(pRest-omega)), 2), 100.);
    if (/*limitMax*T < eLossCut || */limitMax < pow(muqperp_over_T, 2)) return 0.;  

    // Rejection method
    const int max_try = 10000; 
    int n_try = 0; 
    double qperp2_over_T = 0., qperp2_over_T_test, rate_qperp2_test, max_rate = 0., r; 
    // JSINFO << "max rate " << muqperp_over_T; 
    max_rate = Interpolator_dGamma_domega_qperp2(omega_over_T, muqperp_over_T*muqperp_over_T, process); 
    max_rate *= 1.2; 
    while (n_try < max_try)
    {
        qperp2_over_T_test = pow(muqperp_over_T, 2)+ZeroOneDistribution(*GetMt19937Generator())*(limitMax-pow(muqperp_over_T, 2));
        // JSINFO << "qperp2 test " << qperp2_over_T_test; 
        rate_qperp2_test = Interpolator_dGamma_domega_qperp2(omega_over_T, qperp2_over_T_test, process); 
        // JSINFO << "rate " << rate_qperp2_test; 
        r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
        if (rate_qperp2_test > max_rate) JSWARN << "The sampled rate is larger than the maximum rate we assumed in elastic 2d part!! " << rate_qperp2_test << " " << max_rate;
        if (r < rate_qperp2_test)
        {
            qperp2_over_T = qperp2_over_T_test; 
            break; 
        }
        n_try++; 
        if (n_try == max_try) JSWARN << "cannot find a proper qperp2 in elastic part!!!"; 
    }
    // JSINFO << BOLDWHITE << "number of trials is " << n_try; 
    return sqrt(qperp2_over_T) * T; 
}

double Tequila::Energy_Transfer(double pRest, double T, process_type process)
{
    // Rejection method
    double pRest_over_T = pRest / T; 
    const int max_try = 10000; 
    int n_try = 0; 
    double omega_over_T = 0., omega_over_T_test=0., rate_omega_test, max_rate = 0., max_omega_over_T, r; 
    double sample_omega_over_T_max = std::min(Lambda/(2.*T), pRest_over_T - muqperp_over_T); 
    int j = 0; 
    for (size_t i = 0; i < Nw; i++)
    {
        double i_rate = exp(elasticTable[process].y[i])/(elasticTable[process].x[i]*elasticTable[process].x[i]+1.); 
        if (max_rate < i_rate) { max_rate = i_rate; j = i; }
    }
    // max_rate = Interpolator_dGamma_domega(elas_omega_over_T_pos_min, process); 
    max_rate *= 1.5;
    while (n_try < max_try)
    {
        double rn = ZeroOneDistribution(*GetMt19937Generator()); 
        double elas_omega_over_T_pos_span = sample_omega_over_T_max - elas_omega_over_T_pos_min; 
        double elas_omega_over_T_neg_span = elas_omega_over_T_neg_max - elas_omega_over_T_neg_min; 
        if (rn <= elas_omega_over_T_neg_span/(elas_omega_over_T_neg_span+elas_omega_over_T_pos_span))
            omega_over_T_test = elas_omega_over_T_neg_min + rn/elas_omega_over_T_neg_span*(elas_omega_over_T_neg_span+elas_omega_over_T_pos_span)*(elas_omega_over_T_neg_max-elas_omega_over_T_neg_min); 
        else
            omega_over_T_test = elas_omega_over_T_pos_min + (rn-elas_omega_over_T_neg_span/(elas_omega_over_T_neg_span+elas_omega_over_T_pos_span))/elas_omega_over_T_pos_span*(elas_omega_over_T_neg_span+elas_omega_over_T_pos_span)*(elas_omega_over_T_pos_max-elas_omega_over_T_pos_min); 
        // omega_over_T_test = elas_omega_over_T_min+ZeroOneDistribution(*GetMt19937Generator())*(Lambda/(2.*T)-elas_omega_over_T_min);
	if (process == gg || process == gq || process == qg || process == qq || process == qqp || process == qqb) 
            rate_omega_test = Interpolator_dGamma_domega(omega_over_T_test, process) * (1. - omega_over_T_test/pRest_over_T); 
        else
            rate_omega_test = Interpolator_dGamma_domega(omega_over_T_test, process); 

        /*
        double Q = max(mu_min, omega_over_T_test*T); 
        rate_omega_test *= pow(0.1185*log(Q_Z/Lambda_qcd)/log(Q/Lambda_qcd), 2)/alpha_s/alpha_s; 
        max_rate *= pow(0.1185*log(Q_Z/Lambda_qcd)/log(Q/Lambda_qcd), 2)/alpha_s/alpha_s; */
        if (rate_omega_test > max_rate) JSWARN << "The sampled rate is larger than the maximum rate we assumed in elastic 1d part!! "<< omega_over_T_test << " " << elas_omega_over_T_pos_min << " process is " << process; 
        r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
        // JSINFO << max_rate << " " << rate_omega_test; 
        if (r < rate_omega_test)
        {
            omega_over_T = omega_over_T_test; 
            break; 
        }
        n_try++; 
        // if (n_try > 1000) JSINFO << n_try; 
        if (n_try == max_try) JSWARN << "cannot find a proper omega in elastic part!!!"; 
    }
    return omega_over_T * T; 
}

double Tequila::TransverseMomentum_Transfer_Split(double pRest, double omega, double T, process_type process)
{
    // sampled qperp here is real qperp-sqrt(q^2-q_z^2), but not \tilde{q}_\perp in large angle case. 
    double pRest_over_T = pRest/T;
    double omega_over_T = omega/T; 
    // Notice that here we sample qperp^2 but not qperp, so that the rate of qperp^2 is the largest at smallest qperp^2 (which is not the case for qperp)
    // qperp2_over_T here is actually qperp^2/T^2
    double delta_E_over_2Tqperp2 = pRest_over_T / (4.*omega_over_T*(pRest_over_T-omega_over_T));
    // In theory, the upper limit of qperp is 2(p-omega), but we use a smaller limit since the rate is negligible at large qperp 
    double split_qperp2_over_T_max = std::min({10./delta_E_over_2Tqperp2*T*T, pow((pRest_over_T-omega_over_T), 2)}); 
    double split_qperp2_over_T_min = std::min({0.001/delta_E_over_2Tqperp2*T*T, pow((pRest_over_T-omega_over_T), 2)});
    // if (split_qperp2_over_T_max*T < eLossCut) return 0.;  
    // Rejection method
    const int max_try = 10000; 
    int n_try = 0; 
    double qperp2_over_T = 0., qperp2_over_T_test, rate_qperp2_test, max_rate = 0., r; 

    max_rate = splitdGammadxdqperp2Rate(split_qperp2_over_T_min*delta_E_over_2Tqperp2); 
    max_rate *= 1.5; 
    while (n_try < max_try)
    {
        qperp2_over_T_test = split_qperp2_over_T_min+ZeroOneDistribution(*GetMt19937Generator())*(split_qperp2_over_T_max-split_qperp2_over_T_min);
        rate_qperp2_test = splitdGammadxdqperp2Rate(qperp2_over_T_test*delta_E_over_2Tqperp2); 
        r = max_rate * ZeroOneDistribution(*GetMt19937Generator());
        if (rate_qperp2_test > max_rate) JSWARN << "The sampled rate is larger than the maximum rate we assumed in elastic 2d part!! " << rate_qperp2_test << " " << max_rate; 
        if (r < rate_qperp2_test)
        {
            qperp2_over_T = qperp2_over_T_test; 
            break; 
        }
        n_try++; 
        if (n_try == max_try) JSWARN << "cannot find a proper qperp in elastic part!!! "; 
    }
    // qperp here is real qperp, but not \tilde{q}_\perp as sampled in large angle process
    double qperp = sqrt(qperp2_over_T)*T; 
    return qperp; 
}

double Tequila::splitdGammadxdqperp2Rate(double deltaE_over_2T)
{
    double rate_ratio = -1.*log(1.-exp(-1.*deltaE_over_2T)); 
    return rate_ratio; 
}

double Tequila::xSampling(double pRest, double T, process_type process)
{
    double x = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1.-Lambda/pRest); 
    // double x = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1./2-Lambda/(2.*pRest)); 
    double rate_x = splittingRateOmega(pRest, x, T, process); 
    double x_new, rate_x_new, ratio; 

    for (int step = 0; step < Nsteps; step++)
    {
        x_new = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1.-Lambda/pRest); 
        // x_new = Lambda/(2.*pRest) + ZeroOneDistribution(*GetMt19937Generator()) * (1./2-Lambda/(2.*pRest)); 
        rate_x_new = splittingRateOmega(pRest, x_new, T, process); 
        ratio = rate_x_new / rate_x; 
        double rn = ZeroOneDistribution(*GetMt19937Generator()); 
        if (rn < ratio)
        {
            x = x_new; 
            rate_x = rate_x_new; 
        }
    }
    return x; 
}

bool Tequila::is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

bool Tequila::is_exist (const std::string& name)
{
    struct stat buffer;   
    return (stat (name.c_str(), &buffer) == 0); 
}

double Tequila::splittingF(double x, process_type process)
{
        // only the ratio of the differential rate matters, so we do not have correct factor here. 
	double F = 0.; 
	switch(process)
	{
		case gg_split: F = 16.*CA*CA*(3.+(1.-x)/pow(x, 2)+x/pow(1.-x, 2)-x*(1.-x)); 
			 break; 
		case gq_split: F = 4.*(1.+pow(1.-x, 2))/(pow(x, 2)*(1.-x))*(CF*pow(x, 2)+CA*(1.-x))*6./2; 
		 	 break; 
                // case ggqqb_split: F = CF*((1.-x)/x+x/(1.-x))-CA*(x*x+pow(1.-x, 2)); 
                //          break; 
		case qg_split: F = 8.*CF*(1.+pow(1.-x, 2))/(pow(x, 2)*(1.-x))*(CF*pow(x, 2)+CA*(1.-x));
			 break; 
		case qq_split: F = (4.*CF*((1.+pow(1.-x, 2))/pow(x, 2)+(1.+pow(x, 2))/pow(1.-x, 2))+16.*CF*(CF-CA/2)/x/(1.-x))/2/2; 
			 break; 
		case qqp_split: F = 4.*CF*(1.+pow(1.-x, 2))/pow(x, 2)*4./2; 
			 break; 
		case qqb_split: F = 4.*CF*((1.+pow(1.-x, 2))/pow(x, 2)+pow(x, 2)+pow(1.-x, 2))-16.*CF*(CF-CA/2)*pow(1.-x, 2)/x/2; 
			 break; 
                // case qqbp_split: F = CF*(x*x+pow(1.-x, 2));
                //          break;  
		default: JSWARN << "The process determination is wrong! The process is " << process; 
			 break; 
	}
	return F; 
}

double Tequila::splittingRateOmega(double pRest, double x, double T, process_type process)
{
	return pow(g_hard_elas, 4)*pow(T, 2)/(96.*8.*M_PI*pow(pRest, 2)) * splittingF(x, process); 
}

void Tequila::LoadElasticTables()
{
    IntTabulator inttabulator; 
    gsl_interp2d *interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, Nw+2, Nq+1); 
    // Iterate over all the elastic processes
    for (int iProcess = gg; iProcess <= qqb; iProcess++)
    {
        JSINFO << BOLDYELLOW << "process is " << inttabulator.GetProcessString(iProcess); 
        process_type process = static_cast<process_type>(iProcess); 
        tables iTables; 

        if (!is_exist((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+".dat").c_str())) inttabulator.Tabulator_dGamma_domega_qperp2(path_to_tables, process); 
        std::ifstream table2d_in((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+".dat").c_str()); 
        if (is_empty(table2d_in)) inttabulator.Tabulator_dGamma_domega_qperp2(path_to_tables, process); 
        double z; 
        for (size_t iomega = 0; iomega <= Nw+1; iomega++)
        {
	    if (iomega < Nw/5+1.)
    	        iTables.xa[iomega] = -exp(((double)(Nw/5-iomega))*(log(-elas_omega_over_T_neg_min)-log(-elas_omega_over_T_neg_max))/Nw*5.+log(-elas_omega_over_T_neg_max)); 
	    else
		iTables.xa[iomega] = exp(((double)(iomega-Nw/5-1.))*(log(elas_omega_over_T_pos_max)-log(elas_omega_over_T_pos_min))/Nw*5./4+log(elas_omega_over_T_pos_min));
            for (size_t iqperp = 0; iqperp <= Nq; iqperp++)
            {
                iTables.ya[iqperp] = pow(exp(((double)iqperp)*(log(elas_qperp_over_T_max)-log(muqperp_over_T_0))/Nq+log(muqperp_over_T_0)), 2); 
                table2d_in >> z;
                iTables.rate_qperp2[iomega][iqperp] = z;
                gsl_interp2d_set(interp, &za[iProcess*(Nw+2)*(Nq+1)], iomega, iqperp, log(z)); 
            }
        }
        table2d_in.close(); 

        for (size_t iomega = 0; iomega <= Nw+1; iomega++)
        {
            iTables.x[iomega] = iTables.xa[iomega]; 
            // iTables.y[iomega] = 0.; 
            iqperp0 = (int)((log(muqperp_over_T)-log(muqperp_over_T_0))/(log(elas_qperp_over_T_max)-log(muqperp_over_T_0))*Nq); 
            for (size_t iqperp = iqperp0; iqperp < Nq; iqperp++)
            {
                double dy = pow(exp(((double)(iqperp+1))*(log(elas_qperp_over_T_max)-log(muqperp_over_T_0))/Nq+log(muqperp_over_T_0)), 2) - pow(exp(((double)iqperp)*(log(elas_qperp_over_T_max)-log(muqperp_over_T_0))/Nq+log(muqperp_over_T_0)), 2); 
                iTables.y[iomega] += dy * (iTables.rate_qperp2[iomega][iqperp] + iTables.rate_qperp2[iomega][iqperp+1]) / 2; 
            }
        }

        size_t i_total_rate = 0; 	// total rate iteration is different from x and y, because of the gap around 0
        for (size_t iomega = 0; iomega < Nw+1; iomega++)
        {
	        if (abs(iTables.x[iomega]-elas_omega_over_T_neg_max) <= epsilon && iTables.x[iomega] < elas_omega_over_T_pos_min) continue;  
            double dx = iTables.x[iomega+1] - iTables.x[iomega];
            iTables.total_rate[0][i_total_rate] = iTables.x[iomega+1]; 
            iTables.total_rate_w[0][i_total_rate] = iTables.x[iomega+1]; 
            if (iomega == 0)
            {
                iTables.total_rate[1][i_total_rate] = (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2; 
                iTables.total_rate_w[1][i_total_rate] = (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2 * (iTables.x[iomega] + iTables.x[iomega+1]) / 2; 
            }
            else
            {
                iTables.total_rate[1][i_total_rate] = iTables.total_rate[1][i_total_rate-1] + (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2;
                iTables.total_rate_w[1][i_total_rate] = iTables.total_rate_w[1][i_total_rate-1] + (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2 * (iTables.x[iomega] + iTables.x[iomega+1]) / 2;
            }
            i_total_rate++; 
        }
        for (size_t iomega = 0; iomega <= Nw+1; iomega++)
        {
            iTables.y[iomega] = log(iTables.y[iomega]*(iTables.x[iomega]*iTables.x[iomega]+1)); 
        }
        elasticTable.push_back(iTables); 
    }

    for (int iProcess = GqQg; iProcess <= QbqGg; iProcess++)
    {
        process_type process = static_cast<process_type>(iProcess); 
        tables iTables; 
        std::cout << "i process is " << iProcess << "\n"; 
	    if (!is_exist((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+".dat").c_str())) inttabulator.Tabulator_conversion_dGamma_domega_qperp(path_to_tables, process); 
	    std::ifstream table2d_in_2((path_to_tables+"elastic_rate_table"+inttabulator.GetProcessString(process)+".dat").c_str()); 
	    if (is_empty(table2d_in_2)) inttabulator.Tabulator_conversion_dGamma_domega_qperp(path_to_tables, process); 
	    double z; 
    	for (size_t iomega = 0; iomega <= Nw; iomega++)
    	{
    	    iTables.xa[iomega] = ((double)iomega)*(elas_omega_over_T_pos_max-elas_omega_over_T_neg_min)/Nw+elas_omega_over_T_neg_min; 
	        // double qperpMax = 2.*sqrt(kMax*kMax + kMax*iTables.xa[iomega]); 
	        for (size_t iqperp = 0; iqperp <= Nq; iqperp++)
    	    {
    	        iTables.ya[iqperp] = ((double)iqperp)*(elas_qperp_over_T_max-muqperp_over_T_0)/Nq+muqperp_over_T_0; 
		        table2d_in_2 >> z; 
		        iTables.rate_qperp2[iomega][iqperp] = z; 
		        gsl_interp2d_set(interp, &za[iProcess*(Nw+1)*(Nq+1)], iomega, iqperp, log(z)); 
    	    }
    	}
	    table2d_in_2.close(); 

	    for (size_t iomega = 0; iomega <= Nw; iomega++)
    	{
	        iTables.x[iomega] = iTables.xa[iomega]; 
	        iTables.y[iomega] = 0.;
	        // double qperpMax = 2.*sqrt(kMax*kMax + kMax*iTables.x[iomega]); 
	        iqperp0 = (int)((muqperp_over_T-muqperp_over_T_0)/(elas_qperp_over_T_max-muqperp_over_T_0)*Nq); 
	        for (size_t iqperp = iqperp0; iqperp < Nq; iqperp++)
    	    {
	            double dy = (elas_qperp_over_T_max-muqperp_over_T_0)/Nq; 
		        iTables.y[iomega] += dy * (iTables.rate_qperp2[iomega][iqperp] + iTables.rate_qperp2[iomega][iqperp+1]) / 2; 
	        }
	    }

        size_t i_total_rate = 0; 	// total rate iteration is different from x and y, because of the gap around 0
	    for (size_t iomega = 0; iomega < Nw; iomega++)
    	{
	        if (abs(iTables.x[iomega]-elas_omega_over_T_neg_max) <= epsilon && iTables.x[iomega] < elas_omega_over_T_pos_min) continue;  
    	    double dx = iTables.x[iomega+1] - iTables.x[iomega]; 
            iTables.total_rate[0][i_total_rate] = iTables.x[iomega+1]; 
            // iTables.total_rate_w[0][i_total_rate] = iTables.x[iomega+1]; 
            if (iomega == 0)
            {
                iTables.total_rate[1][i_total_rate] = (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2; 
                // iTables.total_rate_w[1][i_total_rate] = (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2 * (iTables.x[iomega] + iTables.x[iomega+1]) / 2; 
            }
            else
            {
                iTables.total_rate[1][i_total_rate] = iTables.total_rate[1][i_total_rate-1] + (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2;
                // iTables.total_rate_w[1][i_total_rate] = iTables.total_rate_w[1][i_total_rate-1] + (iTables.y[iomega+1] + iTables.y[iomega]) * dx / 2 * (iTables.x[iomega] + iTables.x[iomega+1]) / 2;
            }
            i_total_rate++; 
	    }
        // std::cout << "initial conversion rate " << iTables.total_rate[0][i_total_rate-1] << " " << iTables.total_rate[1][i_total_rate-1] << "\n"; 
        // std::cout << "initial conversion rate " << iTables.total_rate[0][int(i_total_rate/10)-1] << " " << iTables.total_rate[1][int(i_total_rate/10)-1] << "\n"; 
	    elasticTable.push_back(iTables); 
    }
    gsl_interp2d_free(interp); 
}

// Use interpolator for total rate to estimate the rate up to Lambda
double Tequila::interpolatorElasticTotalRate(double p_eval, double T, process_type process)
{
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp *interp = gsl_interp_alloc (gsl_interp_linear, Nw);
    gsl_interp_init (interp, elasticTable[process].total_rate[0], elasticTable[process].total_rate[1], Nw); 
    double result; 
    result = gsl_interp_eval (interp, elasticTable[process].total_rate[0], elasticTable[process].total_rate[1], p_eval, acc);  
    gsl_interp_free (interp); 
    gsl_interp_accel_free (acc); 
    return result * T * pow(g_hard_elas, 4); 
}

// The interpolation of integral(w*dw)
double Tequila::interpolatorElasticEnergyLoss(double p_eval, double T, process_type process)
{
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp *interp = gsl_interp_alloc (gsl_interp_linear, Nw);
    gsl_interp_init (interp, elasticTable[process].total_rate_w[0], elasticTable[process].total_rate_w[1], Nw); 
    double result; 
    result = gsl_interp_eval (interp, elasticTable[process].total_rate_w[0], elasticTable[process].total_rate_w[1], p_eval, acc);  
    gsl_interp_free (interp); 
    gsl_interp_accel_free (acc); 
    return result * T * T * pow(g_hard_elas, 4); 
}

double Tequila::Interpolator_dGamma_domega(double omega, process_type process)
{
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_interp *interp = gsl_interp_alloc (gsl_interp_linear, Nw+2);

    gsl_interp_init (interp, elasticTable[process].x, elasticTable[process].y, Nw+2); 
    double result; 
    result = gsl_interp_eval (interp, elasticTable[process].x, elasticTable[process].y, omega, acc);  
    gsl_interp_free (interp); 
    gsl_interp_accel_free (acc); 
    return exp(result)/(omega*omega+1.); 
}

double Tequila::Interpolator_dGamma_domega_qperp2(double omega, double qperp2, process_type process)
{	
    gsl_interp2d *interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, Nw+2, Nq+1);
    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();
    // gsl_interp2d_init(interp, elasticTable[process].xa, elasticTable[process].ya, elasticTable[process].za, Nw+1, Nq+1);
    gsl_interp2d_init(interp, elasticTable[process].xa, elasticTable[process].ya, &za[process*(Nw+2)*(Nq+1)], Nw+2, Nq+1);
    // if (omega == elasticTable[iProcess].xa[iqperp] && qperp == elasticTable[iProcess].ya[iqperp]) return z; 
    double result; 
    // result = gsl_interp2d_eval(interp, elasticTable[process].xa, elasticTable[process].ya, elasticTable[process].za, omega, qperp, xacc, yacc); 
    result = gsl_interp2d_eval(interp, elasticTable[process].xa, elasticTable[process].ya, &za[process*(Nw+2)*(Nq+1)], omega, qperp2, xacc, yacc); 
    gsl_interp2d_free(interp);
    gsl_interp_accel_free(xacc);
    gsl_interp_accel_free(yacc);
    return exp(result); 
}

void Tequila::allocate_memory_for_radiative_rate_table() {
	// Allocate memory for differential and integrated rate
	rate_ggg_p = new double [nb_points_in_p];
	rate_gqq_p = new double [nb_points_in_p];
	rate_qqg_p = new double [nb_points_in_p];
	//maximum_differential_rate = new double [nb_points_in_p];
	differential_rate_ggg_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_gqq_p_omega_qperp = new double ** [nb_points_in_p];
	differential_rate_qqg_p_omega_qperp = new double ** [nb_points_in_p];
	for(int ip=0;ip<nb_points_in_p; ip++) {
		differential_rate_ggg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_gqq_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		differential_rate_qqg_p_omega_qperp[ip]=new double * [nb_points_in_omega];
		for(int iomega=0;iomega<nb_points_in_omega; iomega++) {
			 differential_rate_ggg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_gqq_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
			 differential_rate_qqg_p_omega_qperp[ip][iomega]=new double [nb_points_in_qperp];
		}
	}	
}

//Load rate differential in incoming parton energy p and radiated parton energy omega (not differential in transverse momentum q)
//The table contains dGamma/(dx domega)= gs^4*use_table(), where use_table() is the function from Moore's code that tabulates the collinear rate
void Tequila::load_differential_rate(const double alpha_s, const double alpha_EM, const int Nf, const std::string location_of_collinear_rates)
{
    // Check if the tabulated rate is already available for these parameters
    // The format is "dGamma_dp_domega_NfX_alphasYYY" where "X" is the number of flavours and "YYY" is the value of alpha_s
    std::stringstream filename;
    filename << location_of_collinear_rates << "dGamma_dp_domega_Nf" << Nf << "_alphaEM" << alpha_EM << "_alphaS" << alpha_s;
    //Open file
    std::ifstream rate_file;
    rate_file.open(filename.str().c_str(),std::fstream::in);

    //Struct [from Moore's code] to store rates
    dGammas Moore_rate_arrays;
    Gamma_info Moore_rate_info;
    // Save file with proper name
    sprintf(Moore_rate_info.in_fname, "%s", filename.str().c_str());

    // If not, run Guy's program given parameters
    if (!rate_file.is_open())
    {
        std::cout << "Pre-tabulated collinear rate not available. Tabulating it now...\n";
        //Pre-initialized in data: Nc, Nf, alpha_s, alpha
        Moore_rate_info.Nf=Nf;
        Moore_rate_info.Nc=3;
        Moore_rate_info.alpha_s=alpha_s;
        Moore_rate_info.alpha=alpha_EM;
        build_table(&Moore_rate_info , &Moore_rate_arrays);

        std::cout << "... and writing it into file \"" << filename.str().c_str() << "\" for future use.\n";
        write_table(&Moore_rate_info , &Moore_rate_arrays);

    }
    else
    {
        std::cout << "Collinear rate available for value of alpha_s, alpha_EM and Nf. Reading from file.\n";

        rate_file.close();
        // Read rate file in temporary array differential in p and omega using Moore's functions
        read_table(&Moore_rate_info,&Moore_rate_arrays);
        std::cout << "Collinear rates read.\n";
    }
    const double gs4=alpha_s*alpha_s*(16*M_PI*M_PI);


    // Populate member array "differential_rate_p_omega_qperp[]" with p and omega rate from Moore's code and Gaussian distribution for qperp
    // Some useful parameters for the (temporary) qperp distrib
    // The (temporary) Gaussian envelope used in the q_perp direction (since the q_perp distribution is not currently known)
    const double envelope_sigma_sqr= 4.*M_PI*alpha_s/3.*(3 + nf*0.5); //Debye mass
    for(int ip=0;ip<nb_points_in_p; ip++)
    {
        const double p_over_T_val=get_p_over_T_from_index(ip);
        for(int iomega=0;iomega<nb_points_in_omega; iomega++)
        {
            const double omega_over_T_val=get_omega_over_T_from_index(ip,iomega);
            //Object "Moore_rate_arrays" does not actually contains the rate, but rather the rate stripped of various factors 
            //Moore's function "use_table()" returns the actual rate, after multiplications by the proper factors
            //Using "use_table()" to do this is kind-of sub-optimal, but speed shouldn't be too much of an issue here
            //Need to multiply by g_s^4 since the factor is stripped from the rate in Moore's program
            const double tmp_rate_ggg=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma_ggg , 2 );
            const double tmp_rate_gqq=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma_gqq , 1 );
            const double tmp_rate_qqg=gs4*use_table(p_over_T_val , omega_over_T_val , Moore_rate_arrays.dGamma, 0 );

            for(int iq=0; iq<nb_points_in_qperp;iq++)
            {
                // q_perp envelope function: Gaussian of width sigma \propto m_D?
                // Integrate[ (2 Pi r) Exp[-r^2/sigma^2]/Sqrt[Pi sigma^2]^2, {r, 0, Infinity}]=1
                //jacobian between E dG/dq^3 and dG/(dqperp domega)???
                const double qperp_over_T=get_qperp_over_T_from_index(iq);
                //const double envelope=exp(-qperp_over_T*qperp_over_T/envelope_sigma_sqr)/(M_PI*envelope_sigma_sqr);
                // Hack to make the integral over q_perp unity
                //const double jacobian=2*M_PI*qperp_over_T;
                const double envelope=1./(qperp_over_T_max()); 
                const double rate_ggg=tmp_rate_ggg*envelope;
                const double rate_gqq=tmp_rate_gqq*envelope;
                const double rate_qqg=tmp_rate_qqg*envelope;
                // if (ip ==0 && iomega == 0) std::cout << "qperp_over_T=" << qperp_over_T << " & rate=" << rate << "\n";
                //assign rate
                differential_rate_ggg_p_omega_qperp[ip][iomega][iq]=rate_ggg;
                differential_rate_gqq_p_omega_qperp[ip][iomega][iq]=rate_gqq;
                differential_rate_qqg_p_omega_qperp[ip][iomega][iq]=rate_qqg;
                // differential_rate_p_omega_qperp[ip][iomega][iq]=envelope;
            }
        }
    }

}

//Find position in array corresponding to value of p_over_T and omega_over_T
//Matches Moore's code (it must!)
double Tequila::get_index_from_omega_over_T(const double p_over_T, const double omega_over_T) {

	double b;
	
	if ( omega_over_T < 2 )
	{
		if ( omega_over_T < -1 )
		{
			if ( omega_over_T < -2 )
				b = 60 + 5*omega_over_T;
			else
				b = 70+10*omega_over_T;
		}
		else
		{
			if ( omega_over_T < 1 )
				b = 80 + 20*omega_over_T;
			else
				b = 90 + 10*omega_over_T;
		}
	}
	else if ( omega_over_T < p_over_T-2 )
	{ /* This is that tricky middle ground. */
		b = 190 - 10*log ( 1.000670700260932956l / 
				( 0.0003353501304664781l + (omega_over_T-2) / (p_over_T-4) ) - 1 );
	}
	else
	{
		if ( omega_over_T < p_over_T+1 )
		{
			if ( omega_over_T < p_over_T-1 )
				b = 290 + 10*(omega_over_T-p_over_T);
			else
				b = 300 + 20*(omega_over_T-p_over_T);
		}
		else
		{
			if ( omega_over_T < p_over_T+2 )
				b = 310 + 10*(omega_over_T-p_over_T);
			else
				b = 320 + 5*(omega_over_T-p_over_T);
		}
	}

	return b;

}


//Find the value of p_over_T and omega_over_T for a given position in array
//Matches Moore's code (it must!)
double Tequila::get_omega_over_T_from_index(const int index_p_over_T, const int index_omega_over_T) {

	//Need p_over_T to get omega_over_T
	double p_over_T=get_p_over_T_from_index(index_p_over_T);
	double omega_over_T;

	if ( index_omega_over_T < 50 )        /* spaced by 0.2  from -12 to -2 */
		omega_over_T = -12 + index_omega_over_T * 0.2;
	else if ( index_omega_over_T < 60 )   /* spaced by 0.1  from -2  to -1 */
		omega_over_T = -2 + (index_omega_over_T-50) * 0.1;
	else if ( index_omega_over_T < 100 )  /* spaced by 0.05 from -1  to +1 */
		omega_over_T = -1 + (index_omega_over_T-60) * 0.05;
	else if ( index_omega_over_T < 110 )  /* spaced by 0.1  from +1  to +2 */
		omega_over_T = 1 + (index_omega_over_T-100) * 0.1;
	else if ( index_omega_over_T < 270 )  /* spaced complicated, +2 to p_over_T-2 */
	{
		omega_over_T = 0.1 * (index_omega_over_T-190);
		omega_over_T = 2 + (p_over_T-4) * ( -0.0003353501304664781l
				+ 1.000670700260932956l / (1+exp(-omega_over_T)) );
	}
	else if ( index_omega_over_T < 280 )  /* spaced by 0.1  from p_over_T-2 to p_over_T-1 */
		omega_over_T = p_over_T - 2 + 0.1 * (index_omega_over_T-270);
	else if ( index_omega_over_T < 320 )  /* spaced by 0.05 from p_over_T-1 to p_over_T+1 */
		omega_over_T = p_over_T + 0.05 * (index_omega_over_T - 300);
	else if ( index_omega_over_T < 330 )  /* spaced by 0.1  from p_over_T+1 to p_over_T+2 */
		omega_over_T = p_over_T + 0.1 * (index_omega_over_T - 310);
	else                   /* spaced by 0.2  from p_over_T+2 to p_over_T+12 */
		omega_over_T = p_over_T + 0.2 * (index_omega_over_T - 320);

	return omega_over_T;

}

//Evaluated integrated rate, with cut-off "omega_over_T_cut" on omega, from differential rate stored in member object "differential_rate_p_omega_qperp[][]"
//Also save the maximum value of the rate
void Tequila::evaluate_integrated_rate(double omega_over_T_cut, double *** differential_rate, double * integrated_rate, process_type process) {

	//current omega/T integration not very good if omega_over_T_cut<0.1
	if (omega_over_T_cut<0.1) std::cout << "Warning: omega/T integration is not very good for omega/T cut-off smaller than 0.1\n";
	//loop over all values of "p"
	for(int ip=0;ip<nb_points_in_p; ip++) {

		const double p_over_T=get_p_over_T_from_index(ip);

		// integrate omega/T from -infinity to omega_over_T_cut, and omega_over_T_cut to p/2T
		// the discretization is not uniform in omega/T
		// let's just use the GSL interpolation routine to integrate, since performance is not an issue here
		double integral_omega=0.0;

		// Maximum point in omega/T necessary
		int pOver2T_cut_pos_position_int=ceil(get_index_from_omega_over_T(p_over_T,p_over_T/2.0))+2;
                if (process == qqg)  pOver2T_cut_pos_position_int=ceil(get_index_from_omega_over_T(p_over_T,p_over_T-omega_over_T_cut))+2; 
		//Arrays to store the omega/T rates
		double * omega_rate_array = new double [pOver2T_cut_pos_position_int];
		double * omega_position_array = new double [pOver2T_cut_pos_position_int];
		// Fill an array with the qperp-integrated values for each omega points
		for(int iomega=0;iomega<pOver2T_cut_pos_position_int; iomega++) {

			omega_position_array[iomega]=get_omega_over_T_from_index(ip,iomega);

			// integrate in q_perp
			double integral_qperp=0.0;

			// track the position of various timesteps to compute the weight of each step properly with a minimum number of call to function "get_qperp_over_T_from_index()"
			// trapeze integration not very good for radial coordinate. might require improvements.
			double prev_qperp=0.0;
			double curr_qperp=get_qperp_over_T_from_index(0);
			double next_qperp;
			for(int iq=0; iq<nb_points_in_qperp;iq++) {	
				next_qperp=get_qperp_over_T_from_index(iq+1);
				// get weight
				const double weight_qperp=non_uniform_trapeze_weights(iq,nb_points_in_qperp,prev_qperp,curr_qperp,next_qperp);
				// disable jacobian for now
				// jacobian for radial integration
				//const double jacobian=2*M_PI*curr_qperp;
				const double jacobian=1.0;

				// add contribution from this timestep
				//integral_qperp+=weight_qperp*jacobian*differential_rate_p_omega_qperp[ip][iomega][iq];
				integral_qperp+=weight_qperp*jacobian*differential_rate[ip][iomega][iq];
				//update positions for q_perp
				prev_qperp=curr_qperp;
				curr_qperp=next_qperp;
			}

			omega_rate_array[iomega]=integral_qperp;

		}
		// initialise GSL interpolator
		gsl_interp * interp = gsl_interp_alloc( gsl_interp_akima, pOver2T_cut_pos_position_int);
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_interp_init(interp,  omega_position_array, omega_rate_array, pOver2T_cut_pos_position_int); 

	// integral in omega from -infinity to -omega_over_T_cut
		integral_omega+=gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array, omega_over_T_min(p_over_T), std::max(omega_over_T_min(p_over_T), -1.*omega_over_T_cut), acc);
		// integral in omega from omega_over_T_cut to p_over_T/2
		double omega_over_T_int_max = p_over_T/2; 
		if (process == qqg)
                {
		    omega_over_T_int_max = p_over_T - omega_over_T_cut; 
		    integral_omega += gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array , std::min(omega_over_T_cut, omega_over_T_int_max), omega_over_T_int_max, acc);
                }
                else if (process == ggg)
		    integral_omega += gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array , std::min(omega_over_T_cut, omega_over_T_int_max), omega_over_T_int_max, acc);
                else if (process == gqq)
		    integral_omega += 2. * gsl_interp_eval_integ(interp, omega_position_array, omega_rate_array , std::min(omega_over_T_cut, omega_over_T_int_max), omega_over_T_int_max, acc);
		// free memory
		gsl_interp_free(interp);
		gsl_interp_accel_free(acc);
		delete [] omega_position_array;
		delete [] omega_rate_array;

		//rate_p[ip]=integral_omega;
		integrated_rate[ip]=integral_omega;
	}
}

// weights for the trapezoid rule on a non-uniform grid
// 1/2 Sum[f[getQfromIndex[position]]Which[0==position,(getQfromIndex[position+1]-getQfromIndex[position]),size_of_array-1==position,(getQfromIndex[position]-getQfromIndex[position-1]),True,(getQfromIndex[position+1]-getQfromIndex[position-1])],{position,0,Qnum-1}]
// 1/2 Sum[f[getQfromIndex[position]]Which[0==position,(next_x-curr_x),size_of_array-1==i,(curr_x-prev_x),True,(next_x-prev_x)],{i,0,size_of_array-1}]
double Tequila::non_uniform_trapeze_weights(int position, int size_of_array, double prev_x, double curr_x, double next_x) {

	double weight=0.5;
	
	if (0 == position) {
		weight*=(next_x-curr_x);	
	}
	else if (size_of_array -1 == position) {
		weight*=(curr_x-prev_x);	
	}
	else {
		weight*=(next_x-prev_x);
	}
	return weight;

}

//Returns T^2 dGamma/(dx domega d^2 qperp)
double Tequila::differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T, double *** differential_rate_p_omega_qperp, double T) {


	//tri-linear interpolation
	//somehow I thought this function would be shorter...

	if ((p_over_T<p_over_T_min())||(p_over_T>p_over_T_max())) return 0.0;
	if ((omega_over_T<omega_over_T_min(p_over_T))||(omega_over_T>omega_over_T_max(p_over_T))||(qperp_over_T>qperp_over_T_max())) return 0.0;

	//first, get position in grid of where rate is located
	const double tmp_p=get_index_from_p_over_T(p_over_T);
	const double tmp_omega=get_index_from_omega_over_T(p_over_T,omega_over_T);
	const double tmp_qperp=get_index_from_qperp_over_T(qperp_over_T);

	const int pos_array_p_over_T=floor(tmp_p);
	const int pos_array_omega_over_T=floor(tmp_omega);
	const int pos_array_qperp_over=floor(tmp_qperp);

	//actual positions of the grid points around the desired value
	const double p_over_T_low_val=get_p_over_T_from_index(pos_array_p_over_T);
	const double omega_over_T_low_val=get_omega_over_T_from_index(pos_array_p_over_T,pos_array_omega_over_T);
	const double qperp_over_T_low_val=get_qperp_over_T_from_index(pos_array_qperp_over);
	const double p_over_T_high_val=get_p_over_T_from_index(pos_array_p_over_T+1);
	const double omega_over_T_high_val=get_omega_over_T_from_index(pos_array_p_over_T,pos_array_omega_over_T+1);
	const double qperp_over_T_high_val=get_qperp_over_T_from_index(pos_array_qperp_over+1);

	//value of the rate at the above gridpoints
	const double v000=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T][pos_array_qperp_over];
	const double v001=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T][pos_array_qperp_over+1];
	const double v010=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T+1][pos_array_qperp_over];
	const double v011=differential_rate_p_omega_qperp[pos_array_p_over_T][pos_array_omega_over_T+1][pos_array_qperp_over+1];
	const double v100=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T][pos_array_qperp_over];
	const double v101=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T][pos_array_qperp_over+1];
	const double v110=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T+1][pos_array_qperp_over];
	const double v111=differential_rate_p_omega_qperp[pos_array_p_over_T+1][pos_array_omega_over_T+1][pos_array_qperp_over+1];

	//fraction of each corner to use
	const double frac_p_over_T=(p_over_T-p_over_T_low_val)/(p_over_T_high_val-p_over_T_low_val);
	const double frac_omega_over_T=(omega_over_T-omega_over_T_low_val)/(omega_over_T_high_val-omega_over_T_low_val);
	const double frac_qperp_over_T=(qperp_over_T-qperp_over_T_low_val)/(qperp_over_T_high_val-qperp_over_T_low_val);

	//get the value
	const double v00=v000*(1-frac_qperp_over_T)+v001*frac_qperp_over_T;
	const double v01=v010*(1-frac_qperp_over_T)+v011*frac_qperp_over_T;
	const double v10=v100*(1-frac_qperp_over_T)+v101*frac_qperp_over_T;
	const double v11=v110*(1-frac_qperp_over_T)+v111*frac_qperp_over_T;

	const double v0=v00*(1-frac_omega_over_T)+v01*frac_omega_over_T;
	const double v1=v10*(1-frac_omega_over_T)+v11*frac_omega_over_T;

	double res=v0*(1-frac_p_over_T)+v1*frac_p_over_T;
        double Q = max(omega_over_T*T, mu_min); 
        res *= pow(alpha_QZ*log(Q_Z/Lambda_qcd)/log(Q/Lambda_qcd), 2)/pow(alphas_hard_inel, 2); 
	return res;

}

////Rate
//double ERateColl::rate(const struct ERateParam &rate_params, const Pythia8::Vec4 &p0, const int &id0)
double Tequila::rate_inel(double energy, double temp, double * rate_p)
{
	//energy/temperature ratio of the incoming jet
	//const double energy_over_T=p0.e()/temp;
	const double energy_over_T=energy/temp;

	//if outside of tabulated range, set to 0 (that is, assume untabulated rate is tiny and can be ignored)
	if ((energy_over_T<p_over_T_min())||(energy_over_T>p_over_T_max())) return 0.0;
	//get the real-valued index
	double a = get_index_from_p_over_T(energy_over_T);
	//get the actual index of the array
	int n_p = int(a);
	//get the remainder
	a -= n_p;
	//get rate for energy/T
	double result = (1-a) * rate_p[n_p] + a * rate_p[n_p+1];

	//a factor of temperature is missing from the integrated rate
	result*=temp;


	return result;

}

//Given a value of parton momentum p, sample the rate in omega and qperp
void Tequila::sample_dgamma_dwdq(double p, double T, double *** differential_rate_p_omega_qperp, double &w, double &q, process_type process) {
	//   double lam = pp.lambda(p0) ;
	int ntry = 0  ;
	const int ntry_max = 1000000; 
	
	//helper variables
	//const double qperp_over_T_val_min=0.0;
	//const double qperp_over_T_val_max=qperp_over_T_max();
	const double p_over_T=p/T;
	const double omega_over_T_neg_min=omega_over_T_min(p_over_T);
	const double omega_over_T_neg_max=-muomega_over_T;
	const double omega_over_T_pos_min=muomega_over_T;
	double omega_over_T_pos_max=p_over_T/2.0;
        if (process == qqg)
            omega_over_T_pos_max = std::max(p_over_T - omega_over_T_pos_min, omega_over_T_pos_min); 
	const double omega_over_T_neg_span=(omega_over_T_neg_max-omega_over_T_neg_min);
	const double omega_over_T_pos_span=(omega_over_T_pos_max-omega_over_T_pos_min);

	double max_rate=maximum_rate_p(p_over_T, differential_rate_p_omega_qperp);
        // if (process == gqq) max_rate *= 50.; 

	while (ntry < ntry_max) {
		//      fRateTable.sample_in_x1range(-GSL_DBL_MAX,lam/(2.*T), rng, params, w, q) ;
		//double r = max_rate*rng(params) ;
		double r = max_rate*ZeroOneDistribution(*GetMt19937Generator()) ;
		//      double v = (1. -  T*w/p0.e())/rratio(w) ;
		//const double qx_over_T_tmp=qperp_over_T_val_max*rng(params);
		//const double qy_over_T_tmp=qperp_over_T_val_max*rng(params);
		//const double q_over_T_test=sqrt(qx_over_T_tmp*qx_over_T_tmp+qy_over_T_tmp*qy_over_T_tmp);

		// With[{tmp = (rnd*(xnegspan + xposspan))},
		// If[tmp > xnegspan, xposmin + (tmp - xnegspan), xnegmin + tmp]
		//  ]
		//const double rnd2=rng(params);
		const double rnd2=ZeroOneDistribution(*GetMt19937Generator());
		const double tmp=rnd2*(omega_over_T_neg_span+omega_over_T_pos_span);
		double omega_over_T_test;
		if (tmp < omega_over_T_neg_span) omega_over_T_test=omega_over_T_neg_min+tmp;
		else omega_over_T_test=omega_over_T_pos_min+(tmp-omega_over_T_neg_span);

		// sample q_over_T, remembering it is a radial variable
		// also, q_over_T must be smaller than abs(omega_over_T_test)
		//const double q_over_T_sqr=omega_over_T_test*omega_over_T_test*rng(params);
		//const double q_over_T_test=sqrt(q_over_T_sqr);
		const double q_over_T_test=0.0;  //let's use qperp/T=0 for now

		double v = differential_rate(p_over_T, omega_over_T_test, q_over_T_test, differential_rate_p_omega_qperp, T);

		//safety check: did we pick the maximum rate correctly??
		if (v > max_rate) {
			std::cout << "Function \"maximum_rate_p()\" apparently does not return the maximum of the rate... This is bad. Fix it.\n";
			std::cout << "current guess for maximum=" << max_rate << " & sampled rate at omega/T="<< omega_over_T_test << "&q_perp/T=" << q_over_T_test << " is " << v << "\n";
			exit(1);
			//assert(false);
		}

		if (r < v) {
			w=omega_over_T_test*T;
			q=q_over_T_test*T;
			return ;
		}
		else {
			ntry++ ;
		}
	}
	w=muomega_over_T*T;
	q=0.0;
	std::cout << "*** ERateColl::sample_dgamma_dwdq *** Failed to find a sample "
		"after "<< ntry  << " iterations! Returning with w,q = " << w<<","<< q << std::endl;
}

