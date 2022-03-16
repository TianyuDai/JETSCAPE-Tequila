#ifndef INTTABULATOR_H
#define INTTABULATOR_H

#include "gsl/gsl_integration.h"
#include <fstream>

// qqb is qqbar->qqbar, qqp is qq'->qq', qqbp is qqbar->q'qbar', qqbgg is qqbar->gg
// gq_inel_conv and qg_inel_conv are conversion processes of inelasic part. 

const int CA = 3; 
const int dA = 8; 
const double CF = 4./3.; 
const int dF = 3;
const int Nf = 3;  

const static size_t Nw = 1000; // number of sites in omega grid in tabulator
const static size_t Nq = 1000; 

const double elas_omega_over_T_pos_max = 100.; 
const double elas_omega_over_T_pos_min = 1.e-2; 
const double elas_omega_over_T_neg_max = -1.e-2; 
const double elas_omega_over_T_neg_min = -10.; // should choose to be -2. The rate at min elas w/T is about 10 times smaller than that at max elas w/T.  
const double elas_qperp_over_T_max = 2*sqrt(elas_omega_over_T_pos_max*elas_omega_over_T_pos_max+elas_omega_over_T_pos_max*elas_omega_over_T_neg_min); 
const double muqperp_over_T_0 = 0.08; 

// enum process_type {gg, gq, qg, qq, qqp, qqb, gg_split, gq_split, ggqqb_split, qg_split, qq_split, qqp_split, qqb_split, qqbp_split, ggg, gqq, qqg, none}; 
enum process_type {gg, gq, qg, qq, qqp, qqb, GqQg, QgGq, GgQbq, QbqGg, ggg, gqq, qqg, gg_split, gq_split, qg_split, qq_split, qqp_split, qqb_split, ggqbq_split, qbqgg_split, none};
	
struct f_params
{
	double f_qperp2;
  	double f_p;  
  	double f_k; 
  	double f_q; 
  	double f_phi; 
  	double f_omega; 
  	double f_muperp; 
  	process_type f_process; 
};

class IntTabulator
{
 private: 

	const double ErrAbs = 1.e-9; 
	const double ErrRel = 1.e-3; 
	const int NWorkSpace = 200; 
	const int fKey = 2; 
	const static size_t nProcess = 22; 
	// std::string ProcessStrings[nProcess] = {"gg", "gq", "qg", "qq", "qqp", "qqb", "gg_split", "gq_split", "ggqqb_split", "qg_split", "qq_split", "qqp_split", "qqb_split", "qqbp_split", "ggg", "gqq", "qqg", "none"}; 
        std::string ProcessStrings[nProcess] = {"gg", "gq", "qg", "qq", "qqp", "qqb", "GqQg", "QgGq", "GgQbq", "QbqGg", "ggg", "gqq", "qqg", "gg_split", "gq_split", "qg_split", "qq_split", "qqp_split", "qqb_split", "ggqbq_split", "qbqgg_split", "none"}; 

 public: 
    gsl_integration_workspace *Space_k;
    gsl_integration_workspace *Space_phi ;

  	IntTabulator();
   	virtual ~IntTabulator();
  	
  	std::string GetProcessString( int enumVal ); 
  	
  	double dGamma_domega_qperp2_forTab(double omega, double qperp2, process_type process); 
  	friend double dGamma_domega_qperp2(double qperp2, void *params); 
  	friend double dGamma_domega_qperp2_k(double k, void *params); 
  	friend double dGamma_domega_qperp2_k_phi_gg(double phi, void *params); 
  	friend double dGamma_domega_qperp2_k_phi_gq(double phi, void *params); 
  	friend double dGamma_domega_qperp2_k_phi_qg(double phi, void *params); 
  	friend double dGamma_domega_qperp2_k_phi_qq(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_qqp(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_qqb(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_qqbgg(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_qqbp(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_GqQg(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_QgGq(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_GgQbq(double phi, void *params); 
	friend double dGamma_domega_qperp2_k_phi_QbqGg(double phi, void *params); 

	void Tabulator_dGamma_domega_qperp2(std::string path, process_type process); 
	void Tabulator_omega_dGamma_domega_qperp2(std::string path, process_type process); 
	
}; 

double dGamma_domega_qperp2(double qperp, void *params); 
double dGamma_domega_qperp2_k(double k, void *params); 
double dGamma_domega_qperp2_k_phi_gg(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_gq(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_qg(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_qq(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_qqp(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_qqb(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_qqbgg(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_qqbp(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_GqQg(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_QgGq(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_GgQbq(double phi, void *params); 
double dGamma_domega_qperp2_k_phi_QbqGg(double phi, void *params); 

#endif

