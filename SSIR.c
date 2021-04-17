/* pomp C snippet file: SSIR */
/* Time: 2021-04-17 18:55:06.657 -0400 */
/* Salt: F5C3D9D05BFD17AABD39D608 */

#include <C:/Users/HP/Documents/R/win-library/4.0/pomp/include/pomp.h>
#include <R_ext/Rdynload.h>

 


/* C snippet: 'rinit' */
#define Beta		(__p[__parindex[0]])
#define mu_IR		(__p[__parindex[1]])
#define mu_RS		(__p[__parindex[2]])
#define eta		(__p[__parindex[3]])
#define rho		(__p[__parindex[4]])
#define N		(__p[__parindex[5]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])
#define H		(__x[__stateindex[3]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 S = nearbyint(eta*N);
                      I = 5;
                      R = 0; 
                      H = 5; 
}

#undef Beta
#undef mu_IR
#undef mu_RS
#undef eta
#undef rho
#undef N
#undef S
#undef I
#undef R
#undef H

/* C snippet: 'step.fn' */
#define Beta		(__p[__parindex[0]])
#define mu_IR		(__p[__parindex[1]])
#define mu_RS		(__p[__parindex[2]])
#define eta		(__p[__parindex[3]])
#define rho		(__p[__parindex[4]])
#define N		(__p[__parindex[5]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])
#define H		(__x[__stateindex[3]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 double dN_SI = rbinom(S, 1-exp(-Beta*I/N*dt));
                      double dN_IR = rbinom(I, 1-exp(-mu_IR*dt));
                      double dN_RS = rbinom(R, 1-exp(-mu_RS*dt));
                      S += dN_RS - dN_SI;
                      I += dN_SI - dN_IR;
                      R += dN_IR - dN_RS;
                      H += dN_SI; 
}

#undef Beta
#undef mu_IR
#undef mu_RS
#undef eta
#undef rho
#undef N
#undef S
#undef I
#undef R
#undef H

/* C snippet: 'rmeasure' */
#define Beta		(__p[__parindex[0]])
#define mu_IR		(__p[__parindex[1]])
#define mu_RS		(__p[__parindex[2]])
#define eta		(__p[__parindex[3]])
#define rho		(__p[__parindex[4]])
#define N		(__p[__parindex[5]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])
#define H		(__x[__stateindex[3]])
#define count		(__y[__obsindex[0]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 count = rnbinom(H,rho); 
}

#undef Beta
#undef mu_IR
#undef mu_RS
#undef eta
#undef rho
#undef N
#undef S
#undef I
#undef R
#undef H
#undef count

/* C snippet: 'dmeasure' */
#define Beta		(__p[__parindex[0]])
#define mu_IR		(__p[__parindex[1]])
#define mu_RS		(__p[__parindex[2]])
#define eta		(__p[__parindex[3]])
#define rho		(__p[__parindex[4]])
#define N		(__p[__parindex[5]])
#define S		(__x[__stateindex[0]])
#define I		(__x[__stateindex[1]])
#define R		(__x[__stateindex[2]])
#define H		(__x[__stateindex[3]])
#define count		(__y[__obsindex[0]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 lik = dnbinom(count,H,rho,give_log); 
}

#undef Beta
#undef mu_IR
#undef mu_RS
#undef eta
#undef rho
#undef N
#undef S
#undef I
#undef R
#undef H
#undef count
#undef lik

/* C snippet: 'toEst' */
#define Beta		(__p[__parindex[0]])
#define mu_IR		(__p[__parindex[1]])
#define mu_RS		(__p[__parindex[2]])
#define eta		(__p[__parindex[3]])
#define rho		(__p[__parindex[4]])
#define N		(__p[__parindex[5]])
#define T_Beta		(__pt[__parindex[0]])
#define T_mu_IR		(__pt[__parindex[1]])
#define T_mu_RS		(__pt[__parindex[2]])
#define T_eta		(__pt[__parindex[3]])
#define T_rho		(__pt[__parindex[4]])
#define T_N		(__pt[__parindex[5]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_N = log(N);
	T_Beta = log(Beta);
	T_mu_IR = log(mu_IR);
	T_mu_RS = log(mu_RS);
	T_rho = logit(rho);
	T_eta = logit(eta); 
}

#undef Beta
#undef mu_IR
#undef mu_RS
#undef eta
#undef rho
#undef N
#undef T_Beta
#undef T_mu_IR
#undef T_mu_RS
#undef T_eta
#undef T_rho
#undef T_N

/* C snippet: 'fromEst' */
#define Beta		(__p[__parindex[0]])
#define mu_IR		(__p[__parindex[1]])
#define mu_RS		(__p[__parindex[2]])
#define eta		(__p[__parindex[3]])
#define rho		(__p[__parindex[4]])
#define N		(__p[__parindex[5]])
#define T_Beta		(__pt[__parindex[0]])
#define T_mu_IR		(__pt[__parindex[1]])
#define T_mu_RS		(__pt[__parindex[2]])
#define T_eta		(__pt[__parindex[3]])
#define T_rho		(__pt[__parindex[4]])
#define T_N		(__pt[__parindex[5]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	N = exp(T_N);
	Beta = exp(T_Beta);
	mu_IR = exp(T_mu_IR);
	mu_RS = exp(T_mu_RS);
	rho = expit(T_rho);
	eta = expit(T_eta); 
}

#undef Beta
#undef mu_IR
#undef mu_RS
#undef eta
#undef rho
#undef N
#undef T_Beta
#undef T_mu_IR
#undef T_mu_RS
#undef T_eta
#undef T_rho
#undef T_N

static int __pomp_load_stack = 0;

void __pomp_load_stack_incr (void) {++__pomp_load_stack;}

void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}

void R_init_SSIR (DllInfo *info)
{
R_RegisterCCallable("SSIR", "__pomp_load_stack_incr", (DL_FUNC) __pomp_load_stack_incr);
R_RegisterCCallable("SSIR", "__pomp_load_stack_decr", (DL_FUNC) __pomp_load_stack_decr);
R_RegisterCCallable("SSIR", "__pomp_rinit", (DL_FUNC) __pomp_rinit);
R_RegisterCCallable("SSIR", "__pomp_stepfn", (DL_FUNC) __pomp_stepfn);
R_RegisterCCallable("SSIR", "__pomp_rmeasure", (DL_FUNC) __pomp_rmeasure);
R_RegisterCCallable("SSIR", "__pomp_dmeasure", (DL_FUNC) __pomp_dmeasure);
R_RegisterCCallable("SSIR", "__pomp_to_trans", (DL_FUNC) __pomp_to_trans);
R_RegisterCCallable("SSIR", "__pomp_from_trans", (DL_FUNC) __pomp_from_trans);
}
