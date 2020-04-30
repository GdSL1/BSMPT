/*
 * WallThicknessLib.h
 *
 *  Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller

		This program is free software: you can redistribute it and/or modify
		it under the terms of the GNU General Public License as published by
		the Free Software Foundation, either version 3 of the License, or
		(at your option) any later version.

		This program is distributed in the hope that it will be useful,
		but WITHOUT ANY WARRANTY; without even the implied warranty of
		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
		GNU General Public License for more details.

		You should have received a copy of the GNU General Public License
		along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 */

#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/models/SMparam.h>
#include <complex>

#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>

namespace BSMPT {
namespace ThermalFunctions {

double JbosonIntegrand(const double& x, const double& k, int diff){
	using std::exp;
	using std::log;
	using std::complex;
	complex<double> res(0,0);
	complex<double> kcomplex(k,0), xcomplex(x,0);
	complex<double> TmpSqrt = std::sqrt(kcomplex*kcomplex+xcomplex);
	if(diff == 0){
		res = kcomplex*kcomplex*log(complex<double>(1,0)-exp(-TmpSqrt));
	}
	else if(diff == 1){
		res = kcomplex*kcomplex *exp(-TmpSqrt)/(complex<double>(2,0)*TmpSqrt*(complex<double>(1,0)-exp(-TmpSqrt)));
	}
	else{
		(void) x;
		(void) k;
	}
	return res.real();
}

double JfermionIntegrand(const double& x, const double& k, int diff){
	using std::exp;
	using std::log;
	using std::complex;
	complex<double> res(0,0);
	complex<double> kcomplex(k,0), xcomplex(x,0);
	complex<double> TmpSqrt = std::sqrt(kcomplex*kcomplex+xcomplex);
	if(diff == 0){
		res = kcomplex*kcomplex*log(complex<double>(1,0)+exp(-TmpSqrt));
	}
	else if(diff == 1){
		res = -kcomplex*kcomplex *exp(-TmpSqrt)/(complex<double>(2,0)*TmpSqrt*(complex<double>(1,0)+exp(-TmpSqrt)));
	}
	else{
		(void) x;
		(void) k;
	}
	return res.real();
}


double JbosonNumericalIntegration(const double& x, int diff){
	const size_t workspaceSize = 1000;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (workspaceSize);
	gsl_function F;
	struct GSLTemp{
		double x;
		int diff;
	};
	GSLTemp dummy;
	dummy.diff = diff;
	dummy.x = x;
	F.function = [](double k, void* p){
		struct GSLTemp * params = (struct GSLTemp *) p;
		return JbosonIntegrand(params->x,k,params->diff);
	};
	F.params = &dummy;
	double error, result;
	double upperLimit = 30;
	gsl_integration_qags(&F,0,upperLimit,0,1e-7,workspaceSize,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

double JfermionNumericalIntegration(const double& x, int diff){
	const size_t workspaceSize = 1000;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (workspaceSize);
	gsl_function F;
	struct GSLTemp{
		double x;
		int diff;
	};
	GSLTemp dummy;
	dummy.diff = diff;
	dummy.x = x;
	F.function = [](double k, void* p){
		struct GSLTemp * params = (struct GSLTemp *) p;
		return JfermionIntegrand(params->x,k,params->diff);
	};
	F.params = &dummy;
	double error, result;
	double upperLimit = 30;
	gsl_integration_qags(&F,0,upperLimit,0,1e-7,workspaceSize,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}


double JfermionInterpolatedLow(const double& x, const int& n, int diff){
	if(x == 0 and diff == 0){
		return -7*pow(M_PI,4)/360.0;
	}
	else if(x==0 and diff == 1){
		return pow(M_PI,2)/24;
	}
	using std::pow;
	using std::log;
	double res = 0;
	double cf = 1.5+2*log(4*M_PI)-2*C_euler_gamma-2*log(4);
	if(diff == 0){
		res = -7*pow(M_PI,4)/360.0;
		res +=pow(M_PI,2)/24 *x;
		res += 1/32.0*pow(x,2) * (log(x) - cf);
		double sum = 0;
		for(int l=2;l<=n;l++){
			double Kl = gsl_sf_doublefact(2*l-3)*gsl_sf_zeta(2*l-1)/(gsl_sf_doublefact(2*l)*(l+1))*(pow(2,2*l-1)-1);
			sum += pow(-x/(4*pow(M_PI,2)),l) * Kl ;
		}
		res += -pow(M_PI,2)*x*sum;
	}
	else if(diff == 1){
		res = pow(M_PI,2)/24.0;
		res += x*(-6*cf+3)/96;
		res += x*log(x)/16;
		double sum = 0;
		for(int l=2;l<=n;l++){
			double Kl = gsl_sf_doublefact(2*l-3)*gsl_sf_zeta(2*l-1)/(gsl_sf_doublefact(2*l)*(l+1))*(pow(2,2*l-1)-1);
			sum += -Kl *pow(-x/4.0,l)*(l+1)*pow(M_PI,2-2*l);
		}
		res += sum;
	}



	return res;
}


double JbosonInterpolatedLow(const double& x, const int& n, int diff){
	if(x == 0 and diff == 0){
		return -pow(M_PI,4)/45.0;
	}
	else if(x==0 and diff == 1){
		return pow(M_PI,2)/12.0;
	}
	using std::pow;
	using std::log;
	using std::sqrt;
	double cb = 1.5+2*std::log(4*M_PI)-2*C_euler_gamma;
	double res = 0;
	if(diff == 0){
		res = -pow(M_PI,4)/45.0;
		res += pow(M_PI,2)*x/12.0;
		res += -M_PI*pow(x,1.5)/6;
		res += -pow(x,2)*(log(x)-cb)/32.0;
		double sum = 0;
		for(int l=2;l<=n;l++){
			double Kl = gsl_sf_doublefact(2*l-3)*gsl_sf_zeta(2*l-1)/(gsl_sf_doublefact(2*l)*(l+1));
			sum += pow(-x/(4*pow(M_PI,2)),l) *Kl;
		}
		res += pow(M_PI,2)*x*sum;
	}
	else if(diff == 1){
		res = pow(M_PI,2);
		res += x*(6*cb-3)/96.0;
		res += -x*log(x)/16.0;
		res += -M_PI*sqrt(x)/4.0;
		double sum = 0;
		for(int l=2;l<=n;l++){
			double Kl = gsl_sf_doublefact(2*l-3)*gsl_sf_zeta(2*l-1)/(gsl_sf_doublefact(2*l)*(l+1));
			sum += Kl *pow(-x/4.0,l)*(l+1)*pow(M_PI,2-2*l);
		}
		res += sum;
	}
	return res;
}

double JInterpolatedHigh(const double& x, const int& n, int diff){
	using std::pow;
	using std::exp;
	using std::sqrt;

	double res = 0;
	if(diff == 0){
		double sum = 0;
		for(int l=0;l<=n;l++){
			double Kl = 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l) ;
			sum += Kl * pow(x,-l/2.0);
		}
		res = -exp(-sqrt(x))*sqrt(M_PI/2 * pow(x,1.5))*sum;
	}
	else if(diff == 1){
		double sum = 0;
		for(int l=0;l<=n;l++)
		{
			double Kl = 1/(std::pow(2,l)*gsl_sf_fact(l)) * gsl_sf_gamma(2.5+l)/gsl_sf_gamma(2.5-l);
			sum += Kl *pow(x,(1-l)/2)*(2*l+2*sqrt(x)-3);
		}
		res = exp(-sqrt(x)) *sqrt(2*M_PI)/(8*pow(x,3.0/4.0)) * sum;
	}
	return res;
}

double JfermionInterpolated(const double& x, int diff){
	double res = 0;
	if(x>= C_FermionTheta) {
		res =  -JInterpolatedHigh(x,3,diff);
	}
	else {
		res = -JfermionInterpolatedLow(x,4,diff);
		if(diff == 0) res += -C_FermionShift;
	}
	return res;
}

double JbosonInterpolated(const double& x, int diff){
	double res = 0;
	if(x>= C_BosonTheta){
		res = JInterpolatedHigh(x,3,diff);
	}
	else if(x>= 0){
		res = JbosonInterpolatedLow(x,3,diff);
		if(diff == 0) res += C_BosonShift;
	}
	else if(x<0 and diff == 0){
		res = JbosonInterpolatedNegative(x);
	}
	return res;
}

double JbosonInterpolatedNegative(const double& x){
	if(x>=0) return 0;
	double PotVal = 0;
	double xprev,fprev,xnext,fnext;
	if(-x >= C_NegLine-2) {
		xprev = NegLinearInt[C_NegLine-2][0];
		xnext = NegLinearInt[C_NegLine-1][0];
		fprev = NegLinearInt[C_NegLine-2][1];
		fnext = NegLinearInt[C_NegLine-1][1];
	}
	else{
		size_t pos = (size_t (-x));
		xprev = NegLinearInt[pos][0];
		fprev = NegLinearInt[pos][1];
		xnext = NegLinearInt[pos+1][0];
		fnext = NegLinearInt[pos+1][1];
	}

	PotVal = (fnext-fprev)/(xnext-xprev) * (x-xprev) + fprev;
	PotVal += C_BosonShift;
	return PotVal;
}


}
}
