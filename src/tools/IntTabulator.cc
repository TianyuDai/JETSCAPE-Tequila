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
}

IntTabulator::~IntTabulator()
{
	gsl_integration_workspace_free(Space_phi); 
	gsl_integration_workspace_free(Space_k); 
}

std::string IntTabulator::GetProcessString(int enumVal)
{
  	return ProcessStrings[enumVal];
}

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

void IntTabulator::Tabulator_dGamma_domega_qperp2(std::string path, process_type process)
{
    std::ofstream table2d_out((path+"elastic_rate_table"+GetProcessString(process)+".dat").c_str()); 
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

double dGamma_domega_qperp2_k(double k, void *params)
{	
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
	double q, k, kp, qperp2, omega; 
	qperp2 = p->f_qperp2;
	k = p->f_k; 
	omega = p->f_omega;  
	kp = k + omega; 
	q = sqrt(qperp2 + omega*omega); 
	t = -1.*qperp2; 
        // s here is actually s/2p
	s = (-1.*t/(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	// s = (-1./(2*q*q))*((k+kp)-cos(phi)*sqrt(4*k*kp+t)); 
	u = -1.*s; 
	C = 1./pow(2.*M_PI, 3)/q/2; 
	M2 = (double)(CA*CA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nB(k)*(1+nB(k+omega))); 
	// M2 = (double)(CA*CA)*2.*(pow(s, 2)+pow(u, 2))*(2.*nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_gq(double phi, void *params)
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
	M2 = (double)(dF*CF*CA/dA*6)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qg(double phi, void *params)
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
	M2 = (double)(CF*CA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nB(k)*(1+nB(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qqp(double phi, void *params)
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
	M2 = (double)(dF*CF*CF/dA*4)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qq(double phi, void *params)
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
	M2 = (double)(dF*CF*CF/dA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

double dGamma_domega_qperp2_k_phi_qqb(double phi, void *params)
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
	M2 = (double)(dF*CF*CF/dA)*(pow(s, 2)+pow(u, 2))/pow(t, 2)*(nF(k)*(1-nF(k+omega))); 
	return C*M2; 
}

