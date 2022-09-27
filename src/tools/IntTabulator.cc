#include "IntTabulator.h"
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

IntTabulator::IntTabulator()
{
	Space_phi = gsl_integration_workspace_alloc(NWorkSpace); 
	Space_k = gsl_integration_workspace_alloc(NWorkSpace); 
	Space_w = gsl_integration_workspace_alloc(NWorkSpace); 
	Space_qperp2 = gsl_integration_workspace_alloc(NWorkSpace); 
}

IntTabulator::~IntTabulator()
{
	gsl_integration_workspace_free(Space_phi); 
	gsl_integration_workspace_free(Space_k); 
	gsl_integration_workspace_free(Space_w); 
	gsl_integration_workspace_free(Space_qperp2); 
}

std::string IntTabulator::GetProcessString(int enumVal)
{
  	return ProcessStrings[enumVal];
}

void IntTabulator::Tabulator_split_Gamma(std::string path, process_type process)
{
    std::ofstream table2d_out((path+"elastic_split_rate_table"+GetProcessString(process)+"_runningCoupling.dat").c_str()); 
    double pRest, T; 
    for (size_t i = 0; i <= NT; i++)
    {
        // pRest = exp(((double)i)*(log(elas_pRest_max)-log(elas_pRest_min))/Nq+log(elas_pRest_min)); 
        T = ((double)i)*(elas_T_max-elas_T_min)/NT+elas_T_min; 
        std::cout << "split " << process << " " << T << "\n"; 
	    for (size_t j = 0; j <= Np; j++)
        {
            // std::cout << i << " " << j << "\n";
            pRest = ((double)j)*(sqrt(elas_pRest_max)-sqrt(elas_pRest_min))/Np+sqrt(elas_pRest_min); 
            pRest = pRest * pRest; 
            // std::cout << "pRest " << pRest << "\n";  
            table2d_out << split_Gamma_forTab(pRest, T, process) << " "; 
        }
        table2d_out << "\n"; 
    }
}

double running_coupling(double Q, double T)
{
    double g_running, alpha_s, lambda_QCD=0.2; 
    // double mu_med = 8.*0.3; 
    double Q_med = 4.*T;
    double MZ = 91.1876, alpha_s0 = 0.1189; 
    double c = alpha_s0 / (4.*M_PI/9/std::log(std::pow(MZ, 2)/pow(lambda_QCD, 2)));  
    alpha_s = 4.*M_PI/9/std::log(std::pow(std::max(Q*T, Q_med), 2)/pow(lambda_QCD, 2)) * c;
    // alpha_s = 0.3; 
    g_running = std::sqrt(4.*M_PI*alpha_s); 
    return g_running;  
}

double IntTabulator::split_Gamma_forTab(double pRest, double T, process_type process)
{
	struct f_params p; 
	p.f_pRest = pRest;
    p.f_T = T;  
	p.f_process = process; 
	gsl_function F; 
	F.function = splittingRateOmega; 
	F.params = &p; 
	double result, err; 
    double pcut = 2.; 
    double Lambda = std::min({pcut, std::sqrt(3.*pRest*T), pRest/2});
    // std::cout << "pRest and T " << pRest << " " << T << "\n"; 
    if (pRest/2 > Lambda)
    {
        if (process == gg_split || process == qqb_split) 
            gsl_integration_qag(&F, Lambda, pRest/2, ErrAbs, ErrRel, NWorkSpace, fKey, Space_w, &result, &err); 
        else
            gsl_integration_qag(&F, Lambda, pRest-Lambda, ErrAbs, ErrRel, NWorkSpace, fKey, Space_w, &result, &err); 
        // std::cout << "Lambda and pRest " << Lambda << " " << pRest << " " << result << "\n"; 
        return result; 
    }
    return 0; 
}


double splittingRateOmega(double omega, void *params)
{
    // std::cout << "omega " << omega << "\n"; 
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	// struct f_params p; 
	// IntTabulator cls; 
    double pRest, T; 
    process_type process; 
    pRest = p->f_pRest; 
    T = p->f_T; 
    process = p->f_process; 
	double dGdw = pow(running_coupling(omega/T, T), 4)*T*T/(96.*8.*M_PI*pow(pRest, 2)) * cls.splittingF(omega/pRest, process);
    // std::cout << "T " << T << " " << cls.splittingF(omega/pRest, process) << " " << dGdw << "\n";  
    // if (process == gg_split || process == qqb_split) 
    //     dGdw *= 2.; 
    return dGdw; 
	// return pow(g_hard_elas, 4)*pow(T, 2)/(96.*8.*M_PI*pow(pRest, 2)) * splittingF(x, process); 
}


double IntTabulator::splittingF(double x, process_type process)
{
    // only the ratio of the differential rate matters, so we do not have correct factor here. 
	double F = 0.; 
	switch(process)
	{
		case gg_split: F = 16.*CA*CA*(3.+(1.-x)/pow(x, 2)+x/pow(1.-x, 2)-x*(1.-x)); 
			 break; 
		case gq_split: F = 4.*(1.+pow(1.-x, 2))/(pow(x, 2)*(1.-x))*(CF*pow(x, 2)+CA*(1.-x))*6./2; 
		 	 break; 
        case ggqqb_split: F = CF*((1.-x)/x+x/(1.-x))-CA*(x*x+pow(1.-x, 2)); 
             break; 
		case qg_split: F = 8.*CF*(1.+pow(1.-x, 2))/(pow(x, 2)*(1.-x))*(CF*pow(x, 2)+CA*(1.-x));
			 break; 
		case qq_split: F = (4.*CF*((1.+pow(1.-x, 2))/pow(x, 2)+(1.+pow(x, 2))/pow(1.-x, 2))+16.*CF*(CF-CA/2)/x/(1.-x))/2/2; 
			 break; 
		case qqp_split: F = 4.*CF*(1.+pow(1.-x, 2))/pow(x, 2)*4./2; 
			 break; 
		case qqb_split: F = 4.*CF*((1.+pow(1.-x, 2))/pow(x, 2)+pow(x, 2)+pow(1.-x, 2))-16.*CF*(CF-CA/2)*pow(1.-x, 2)/x/2; 
			 break; 
        case qqbp_split: F = CF*(x*x+pow(1.-x, 2));
             break;  
		default: std::cout << "The process determination is wrong! The process is " << process << "\n"; 
			 break; 
	}
	return F; 
}

void IntTabulator::Tabulator_dGamma_domega(std::string path, process_type process)
{
    std::ofstream table2d_out((path+"elastic_rate_table"+GetProcessString(process)+"_runningCoupling.dat").c_str()); 
    double wp, w, T; 
    for (size_t j = 0; j <= NT; j++)
    {
        T = (elas_T_max - elas_T_min) / NT * j + elas_T_min; 
        for (size_t i = 0; i <= Nw+1; i++)
        {
            if (i < Nw/5+1.)
            {
                wp = (Nw/5-(double)i)*(log(-elas_omega_over_T_neg_min)-log(-elas_omega_over_T_neg_max))/Nw*5.+log(-elas_omega_over_T_neg_max); 
                w = -exp(wp); 
            }
            else
            {
                wp = ((double)(i-Nw/5-1.))*(log(elas_omega_over_T_pos_max)-log(elas_omega_over_T_pos_min))/Nw*5./4+log(elas_omega_over_T_pos_min); 
                w = exp(wp); 
            }
            table2d_out << dGamma_domega_forTab(w, T, process) << " "; 
            // std::cout << "T " << T << " w " << w << " " << dGamma_domega_forTab(w, T, process) << "\n";  
        }
        std::cout << "T " << T << " w " << w << " " << dGamma_domega_forTab(w, T, process) << "\n";  
        table2d_out << "\n"; 
    }
}

double IntTabulator::dGamma_domega_forTab(double omega, double T, process_type process)
{
    struct f_params p;
    p.f_omega = omega;
    p.f_T = T;  
    p.f_process = process; 
    gsl_function F; 
    F.function = dGamma_domega_qperp2; 
    F.params = &p; 
    double result, err;
    // std::cout << "omega and T " << omega << " " << T << "\n";  
	gsl_integration_qagiu(&F, muqperp_over_T*muqperp_over_T, ErrAbs, ErrRel, NWorkSpace, Space_qperp2, &result, &err); 
	// gsl_integration_qagiu(&F, muqperp_over_T, ErrAbs, ErrRel, NWorkSpace, Space_qperp2, &result, &err); 
	// gsl_integration_qag(&F, muqperp_over_T*muqperp_over_T, 125., ErrAbs, ErrRel, NWorkSpace, fKey, Space_qperp2, &result, &err); 
	return result; 
    
}

/*
double IntTabulator::dGamma_domega_qperp2_forTab(double omega, double qperp2, process_type process)
{
	struct f_params p; 
	p.f_qperp2 = qperp2; 
	p.f_omega = omega; 
	p.f_process = process; 
	gsl_function F; 
	F.function = dGamma_domega_qperp2_k; 
	F.params = &p; 
	double result, err; 
	double q = sqrt(qperp2 + omega*omega); 
	double lowLimit = (q - omega) / 2.; 
	gsl_integration_qagiu(&F, lowLimit, ErrAbs, ErrRel, NWorkSpace, Space_k, &result, &err); 
	return result; 
}

// for conversion process only? 
double IntTabulator::dGamma_domega_qperp_forTab(double omega, double qperp, process_type process)
{
	struct f_params p; 
	p.f_qperp2 = qperp; 
	p.f_omega = omega; 
	p.f_process = process; 
	gsl_function F; 
	F.function = dGamma_domega_qperp2_k; 
	F.params = &p; 
	double result, err; 
	double q = sqrt(qperp*qperp + omega*omega); 
	double lowLimit = (q - omega) / 2.; 
	gsl_integration_qagiu(&F, lowLimit, ErrAbs, ErrRel, NWorkSpace, Space_k, &result, &err); 
	return result; 
}

void IntTabulator::Tabulator_dGamma_domega_qperp2(std::string path, process_type process)
{
    std::ofstream table2d_out((path+"elastic_rate_table"+GetProcessString(process)+"_runningCoupling.dat").c_str()); 
    double wp, w, q, qp; 
    for (size_t i = 0; i <= Nw+1; i++)
    {
        if (i < Nw/5+1.)
        {
            wp = (Nw/5-(double)i)*(log(-elas_omega_over_T_neg_min)-log(-elas_omega_over_T_neg_max))/Nw*5.+log(-elas_omega_over_T_neg_max); 
            w = -exp(wp); 
        }
        else
        {
            wp = ((double)(i-Nw/5-1.))*(log(elas_omega_over_T_pos_max)-log(elas_omega_over_T_pos_min))/Nw*5./4+log(elas_omega_over_T_pos_min); 
            w = exp(wp); 
        }

	    for (size_t j = 0; j <= Nq; j++)
        {
            qp = ((double)j)*(log(elas_qperp_over_T_max)-log(muqperp_over_T_0))/Nq+log(muqperp_over_T_0); 
            q = exp(qp);
            double q2 = q*q; 
            table2d_out << dGamma_domega_qperp2_forTab(w, q2, process) << " "; 
        }
        table2d_out << "\n"; 
    }
}
*/
void IntTabulator::Tabulator_conversion_dGamma_domega(std::string path, process_type process)
{
    std::ofstream table2d_out((path+"elastic_rate_table"+GetProcessString(process)+"_runningCoupling.dat").c_str()); 
    double wp, w, T; 
    // flat tabulation (not log) is enough for integration of wdGdw
    for (size_t j = 0; j < NT; j++)
    {
        T = (elas_T_max - elas_T_min) * j + elas_T_min; 
        for (size_t i = 0; i <= Nw; i++)
        {
            w = ((double)i)*(elas_omega_over_T_pos_max-elas_omega_over_T_neg_min)/Nw+elas_omega_over_T_neg_min; 
	        table2d_out << dGamma_domega_forTab(w, T, process) << " "; 
	    }
        table2d_out << "\n"; 
    }
}

double dGamma_domega_qperp2(double qperp2, void *params)
{	
    // std::cout << "qperp2 and lowLimit " << qperp2 << "\n";  
	struct f_params *p = (struct f_params *)params;
    double omega = p->f_omega;  
	IntTabulator cls;
    // qperp2 = 4.;  
	p->f_qperp2 = qperp2; 
	gsl_function F; 
    F.function = dGamma_domega_qperp2_k; 
	F.params = p; 
	double result, err; 
	double q = sqrt(qperp2 + omega*omega); 
	// double q = sqrt(qperp2*qperp2 + omega*omega); 
	double lowLimit = (q - omega) / 2.;
    // std::cout << "omega, qperp2 and lowLimit " << omega << " " << qperp2 << " " << lowLimit << "\n";  
	gsl_integration_qagiu(&F, lowLimit, cls.ErrAbs, cls.ErrRel, cls.NWorkSpace, cls.Space_k, &result, &err);
    // std::cout << omega << " " << qperp2 << " " << result << "\n";  
	return result; 
}

double dGamma_domega_qperp2_k(double k, void *params)
{
    // std::cout << "k start\n"; 
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	p->f_k = k; 
	gsl_function F; 
	switch(p->f_process)
	{
		case gg: F.function = dGamma_domega_qperp2_k_phi_gg; 
			 break; 
		case gq: F.function = dGamma_domega_qperp2_k_phi_gq; 
			 break; 
		case qg: F.function = dGamma_domega_qperp2_k_phi_qg; 
			 break; 
		case qq: F.function = dGamma_domega_qperp2_k_phi_qq; 
			 break; 
		case qqp: F.function = dGamma_domega_qperp2_k_phi_qqp; 
			 break; 
		case qqb: F.function = dGamma_domega_qperp2_k_phi_qqb; 
			 break; 
		case GqQg: F.function = dGamma_domega_qperp_k_phi_GqQg; 
			 break; 
		case QgGq: F.function = dGamma_domega_qperp_k_phi_QgGq; 
			 break; 
		case GgQbq: F.function = dGamma_domega_qperp_k_phi_GgQbq; 
			 break; 
		case QbqGg: F.function = dGamma_domega_qperp_k_phi_QbqGg; 
			 break; 
		default: std::cout << "The process determination is wrong! The process is " << p->f_process; 
			 break; 
	}
	F.params = p; 
	double result, err; 
	gsl_integration_qag(&F, 0, 2.*M_PI, cls.ErrAbs, cls.ErrRel, cls.NWorkSpace, cls.fKey, cls.Space_phi, &result, &err); 
	return result/(2.*M_PI); 
}

double dGamma_domega_qperp2_k_phi_gg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2;
	k = p->f_k; 
	omega = p->f_omega;
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	// q = sqrt(qperp2*qperp2 + omega*omega); 
	// t = -1.*qperp2*qperp2; 
        // s here is actually s/2p
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	// s = (-1./(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
	// C = 1./pow(2.*M_PI, 3)/q*qperp2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4);
    // std::cout << "running coupling " << std::pow(running_coupling(std::sqrt(qperp2), T), 4) << "\n"; 
	M2 = (double)(CA*CA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nB(k)*(1+nB(k+omega)));
	// M2 = (double)(CA*CA)*2.*(pow(s, 2)+pow(u, 2))*(2.*nB(k)*(1+nB(k+omega)));
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_gq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega; 
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(dF*CF*CA/dA*6)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(CF*CA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qqp(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(dF*CF*CF/dA*4)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(dF*CF*CF/dA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qqb(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(dF*CF*CF/dA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_GqQg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	// C = 1./pow(2.*M_PI, 3)/q*qperp; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
        // For conversion process, we neglect T/p here. The additional factor of 1/2 comes from 2p/4p^2. There should be a factor of 2p in s, so for large-angle interactions, the term is s^2, which makes the factor be canceled. 
	M2 = (double)(CF*CF)*u/t*(nF(k)*(1+nB(k+omega)))/2.*2.*Nf;
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_QgGq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(CF*CF)*u/t*(nB(k)*(1-nF(k+omega)))/2.; 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_GgQbq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(CF*CF)*u/t*(nB(k)*(1-nF(k+omega)))/2.*2.*Nf; 
	return C*M2; 
}

double dGamma_domega_qperp_k_phi_QbqGg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega, T; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
    T = p->f_T; 
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
    C *= std::pow(running_coupling(std::sqrt(qperp2), T), 4); 
	M2 = (double)(CF*CF)*u/t*(nF(k)*(1+nB(k+omega)))/2.; 
	return C*M2; 
}
/*
double dGamma_domega_qperp2_k_phi_GqQg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
        // For conversion process, we neglect T/p here. The additional factor of 1/2 comes from 2p/4p^2. There should be a factor of 2p in s, so for large-angle interactions, the term is s^2, which makes the factor be canceled. 
	M2 = (double)(CF*CF)*u/t*(nF(k)*(1+nB(k+omega)))/2.*2.*Nf; 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_QgGq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
	M2 = (double)(CF*CF)*u/t*(nB(k)*(1-nF(k+omega)))/2.; 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_GgQbq(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
	M2 = (double)(CF*CF)*u/t*(nB(k)*(1-nF(k+omega)))/2.*2.*Nf; 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_QbqGg(double phi, void *params)
{
	struct f_params *p = (struct f_params *)params; 
	IntTabulator cls; 
	double s, t, u, M2, C; 
	double q, k, kp, qperp2, omega; 
	qperp2 = p->f_qperp2; 
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
	M2 = (double)(CF*CF)*u/t*(nF(k)*(1+nB(k+omega)))/2.; 
	return C*M2; 
}
*/
