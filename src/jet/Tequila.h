#ifndef TEQUILA_H
#define TEQUILA_H

#include "JetEnergyLossModule.h"
#include "JetScapeConstants.h"
#include "gsl/gsl_integration.h"
#include "IntTabulator.h"
#include "lorentz.h"
#include <vector>
#include "random.h"
#include <math.h>
#include <gsl//gsl_spline.h>
#include "Tequila/evolve9_text_noExp.h"

using namespace Jetscape;

class Tequila : public JetEnergyLossModule<Tequila> //, public std::enable_shared_from_this<Tequila>
{
    static Pythia8::Pythia InternalHelperPythia;
    private: 
    // double alphas_soft;
    // double alphas_hard_elas;
    double beta_para, beta_perp; 
    double alphas_hard_inel;
	// double g_soft;
	// double g_hard_elas;
	double g_hard_inel;
    double T_star; 
   	double M = 0.; 
   	double muqperp_over_T;
    const double eLossCut = 2.;  
   	// AMY rates are calculated in p/T > AMY_p_over_T_cut
  	static constexpr double AMY_p_over_T_cut = 4.01;
    static int pLabelNew;

  	double Q0;
  	double alpha_EM;
  	double hydro_Tc;
  	double muomega_over_T;
    // double qhat_coef; 
	double Lambda; 
    int recoil_on;
    double tStart = 0.1; 
    double pcut;  

  	const double nc = 3.; 
	const double ln2 = log(2); 
    const double epsilon = 1e-8;
    const double mu_min = 2.; 
    const double Lambda_QCD = 0.2; 
    const double alpha_QZ = 0.1185; 
    const double Q_Z = 90.2; 
         

	const int Nsteps = 500; 
	const static size_t nElasProcess = 6;
    const static size_t nConvProcess = 0;
    const static size_t nSplitProcess = 8;
    const static size_t nInelProcess = 3;
    const static size_t nTotalProcess = nElasProcess + nConvProcess + nSplitProcess + nInelProcess + 1; 
	// const static size_t nTotalProcess = 18; 
	
	int iqperp0 = 0; 
	double max_omega_rate; 
	double max_qperp2_rate; 
	double casimir[nElasProcess] = {CA*CA/(24.*M_PI), dF*CF*CA/dA/(8.*M_PI), CF*CA/(24.*M_PI), dF*CF*CF/dA/(72.*M_PI), dF*CF*CF/dA/(18.*M_PI), dF*CF*CF/dA/(72.*M_PI)};		// The coefficient for total rate integration of large omega

	// gg, gq, qg, qq, qqp, qqb
	// There are two more processes in splitting approximation than large-angle elastic interaction
/*
        double splitting_c_lambda[nSplitProcess] = {4.*CA*CA, 16.*CA*2.*nf, 4.*CF*CA, 16.*CF, 16.*CF*(2.*nf-2.), 16.*CF, 0., 0.};
        double splitting_c_p[nSplitProcess] = {5./6*CA*CA, 2.*nf*(2.*CF-4.*CA), CF*CF/2-CF*CA, -4.*CF, -4.*CF*(2.*nf-2), 4.*CF*(6.*CF-3.*CA-1./3), nf*(-CA/3-CF), -8.*CF/3*(CA+3.*CF)};
        double splitting_c_ln[nSplitProcess] = {-2.*CA*CA, 2.*nf*(4.*CF-8.*CA), CF*CF-2.*CF*CA, 8.*CF*(2.*CF-CA-1.), -8.*CF*(2.*nf-2.), -8.*CF*(2.*CF-CA+1.), nf*CF, 8.*CF};
        double Toverp_c_ln[nElasProcess+nConvProcess] = {-2.*CA*CA, -CA*CA, -2.*CF*CA, CF*(2.*CF-CA-1.), -CF, -CF*(2.*CF-CA+1.), -2.*nf*CF*CF/2, -CF*CF, -2.*nf*CF*CF, -CF*CF/2};
*/
	struct tables
	{

		// double total_rate[2][Nw+1]; 
		// double total_rate_w[2][Nw+1]; 
		double dGdw[NT+1][Nw+2]={0.}, dGdq[Nw+1]={0.}; 
		double w[Nw+2], w_tmp[Nw+1]; 
		double T_tmp[NT+1]; 
        double total_rate[NT+1][Nw+1];  
	}; 
	vector <tables> elasticTable; 
	double *za = (double*) malloc((nElasProcess+nConvProcess) * (NT+1) * (Nw+2) * sizeof(double)); 
	double *total_rate_save = (double*) malloc((nElasProcess+nConvProcess) * (NT+1) * (Nw+1) * sizeof(double)); 
	// double *dGdw_save = (double*) malloc((nElasProcess+nConvProcess) * (NT+1) * (Nw+2) * sizeof(double)); 

	struct split_tables
	{
		// double total_rate[2][Nw+1]; 
		// double dGdw[Nw+1]={0.}, dGdq[Nw+1]={0.}; 
		// double w[Nw+2], q[Nq+1]; 
		double p_tmp[Np+1], T_tmp[NT+1]; 
		// double dGdw[Np+1]; 
        double total_rate[NT+1][Np+1];  
	}; 
	vector <split_tables> splitTable; 
	double *zb = (double*) malloc((nSplitProcess) * (NT+1) * (Np+1) * sizeof(double)); 
	double *split_rate_save = (double*) malloc((nSplitProcess) * (NT+1) * (Np+1) * sizeof(double)); 
    
	struct fixed_tables
	{
		// double x[Nw_fixed+2], y[Nw_fixed+1]={0.}; 
		double w[Nw_fixed+2], q[Nq+1]; 
		double rate_qperp2[Nw_fixed+2][Nq+1];  
	}; 
	vector <fixed_tables> elastic_fixedCoup_Table; 
	double *zc = (double*) malloc((nElasProcess+nConvProcess) * (Nw_fixed+2) * (Nq+1) * sizeof(double)); 
    
	bool is_empty(std::ifstream& pFile); 
	bool is_exist (const std::string& name); 
	
	// size of the array containing the differential rate in p, omega and q_\perp
  	// those number in p and omega are, and must be, the same used in Guy Moore's code that is used to evaluate the collinear rates
	// the number of point in q_\perp is arbitrary at this point
	static const int nb_points_in_p=NP, nb_points_in_omega=NK;
	// Addionnal grid in q_perp, not defined in Guy Moore's code
	static const int nb_points_in_qperp=4;
	// the minimum and maximum parton energy p and parton energy loss omega are set in Guy Moore's code as well and must be the same here
	double p_over_T_min() { return 4.01; };
	double p_over_T_max() { return (p_over_T_min()*pow(1.04119,nb_points_in_p-1)); };
	double omega_over_T_min(const double p_over_T) { return -12.0; };
	double omega_over_T_max(const double p_over_T) { return p_over_T + 0.2 * (nb_points_in_omega -1 - 320); };
	// for future use
	double qperp_over_T_max() { return 10.0; };
	//given the above tabulation of p/T, this function returns the index for a given p/T
	//the index is a real number, such that "int(index)" gives the position in the array and "index-int(index)" gives the remainder
	double get_index_from_p_over_T(const double p_over_T) { return 24.7743737154026 * log( p_over_T * .2493765586034912718l ); };
	double get_index_from_omega_over_T(const double p_over_T, const double omega_over_T); //Big function, defined elsewhere instead
	double get_index_from_qperp_over_T(const double qperp_over_T)
        { 
            const double qmin=0.0;
            const double qmax=qperp_over_T_max();
            return (((nb_points_in_qperp-1)*(qperp_over_T-qmin))/(qmax-qmin));
	}
	double get_p_over_T_from_index(const int index_p_over_T) { return (p_over_T_min()*pow(1.04119,index_p_over_T)); };
	double get_omega_over_T_from_index(const int index_p_over_T, const int index_omega_over_T); //Big function, defined elsewhere instead
	// Assume uniform grid in q_perp for now
	double get_qperp_over_T_from_index(const int index_qperp_over_T)
        { 
            const double qmin=0.0;
            const double qmax=qperp_over_T_max();
            return qmin+(qmax-qmin)/(nb_points_in_qperp-1)*index_qperp_over_T;
	};
	// arrays containing the differentail and integrated rates
	double *** differential_rate_qqg_p_omega_qperp, *** differential_rate_ggg_p_omega_qperp, *** differential_rate_gqq_p_omega_qperp;
	double *rate_qqg_p, * rate_ggg_p, *rate_gqq_p;
	//maximum of the differential rate, used to sample the rate

	void allocate_memory_for_radiative_rate_table();

	double non_uniform_trapeze_weights(int position, int size_of_array, double prev_x, double curr_x, double next_x);

	// load differential rate from file into member array "differential_rate_p_omega_q"
	void load_differential_rate(const double alpha_s, const double alpha_EM, const int Nf, const std::string location_of_collinear_rates);
	// fill member array "rate_p" by integrating the content of "differential_rate_p_omega_q"
	void evaluate_integrated_rate(double omega_over_T_cutoff, double *** differential_rate_gqq_p_omega_qperp, double * rate_gqq_p, process_type process);

	double differential_rate(const double p_over_T, const double omega_over_T, const double qperp_over_T, double *** differential_rate_p_omega_qperp, double T);
	double rate_inel(double energy, double temp, double * rate_p);
	void sample_dgamma_dwdq(double p, double T, double *** differential_rate_p_omega_qperp, double &w, double &q, process_type process);
	//educated guess on the highest value of the rate
	//double *maximum_differential_rate;
	double maximum_rate_p(double p_over_T, double *** differential_rate_p_omega_qperp)
        {
            //the collinear rate should presumably be maximum at qperp/T=0 and either omega/T=+/-omega_over_T_cut 
            //assume that the temperature is 0.3, might have some problem if the temperature is too different... 
            const double val1=2.*differential_rate(p_over_T,muomega_over_T,0,differential_rate_p_omega_qperp, 0.3);
            const double val2=2.*differential_rate(p_over_T,-1*muomega_over_T,0,differential_rate_p_omega_qperp, 0.3);
            return val1 > val2 ? val1 : val2;
	};

    // Allows the registration of the module so that it is available to be used by the Jetscape framework.
    static RegisterJetScapeModule<Tequila> reg;

    public:

	Tequila();
	virtual ~Tequila();

  	void LoadElasticTables(); 
    double interpolatorElasticTotalRate(double p_eval, double T_eval, process_type process); 
    double interpolatorElasticTotalRate2d(double p_eval, double T_eval, process_type process); 
    double interpolatorSplitTotalRate(double p_eval, double T_eval, process_type process); 
    double interpolatorSplitTotalRate2d(double p_eval, double T_eval, process_type process); 
    // double interpolatorElasticEnergyLoss(double p_eval, double T, process_type process); 
	double Interpolator_dGamma_domega(double omega, double T_eval, process_type process); 
	// double Interpolator_dGamma_domega_2d(double omega, double T_eval, process_type process); 
	double Extrapolator_dGamma_domega(double omega, process_type process); 
	double Interpolator_dGamma_domega_qperp2(double omega, double qperp2, process_type process); 
	double splittingF(double x, process_type process); 
	double splittingRateOmega(double pRest, double x, double T, process_type process);
    double running_coupling(double Q, double T);  
	// double splittingRateQperp(double pRest, double qperp, double omega, double T, process_type process); 

  	void Init();
	void Finish(); 
  	void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  	process_type DetermineProcess(double pRest, double T, double deltaTRest, int Id); 
  	FourVector Momentum_Update(double omega, double qperp, double T, FourVector pVec); 
  	double TransverseMomentum_Transfer(double pRest, double omega, double T, process_type process); 
  	double Energy_Transfer(double pRest, double T, process_type process); 
	double TransverseMomentum_Transfer_Split(double pRest, double omega, double T, process_type process);
    double splitdGammadxdqperp2Rate(double deltaE_over_2T);  
  	double Energy_Transfer_Split(double pRest, double T, process_type process); 
	double xSampling(double pRest, double T, process_type process); 
	double thermalDistributionSampling(process_type process);
    FourVector getThermalVec(FourVector qVec, double T, int id); 
    double getThermal(double k_min, double T, int kind); 
  	
  	double qhatpara(double E, double T, int id); 
  	double qhatperp(double E, double T, int id); 
  	FourVector Langevin_Update(double dt, double T, FourVector pIn, int Id); 

    static void IncrementpLabel() {pLabelNew++;}
  	
    protected:
  	uniform_real_distribution<double> ZeroOneDistribution;
	string path_to_tables;
}; 

double dGamma_domega(double omega, void *params); 
double dGamma_domega_qperp2(double qperp2, void *params); 
double dGamma_domega_qperp2_k(double k, void *params); 
double dGamma_domega_qperp2_k_phi(double phi, void *params); 

#endif

