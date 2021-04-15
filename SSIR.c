/* pomp C snippet file: SSIR */
/* Time: 2021-04-15 15:28:57.116 -0400 */
/* Salt: D30BE1DD3CB1B38C56FEA108 */

#include <C:/Users/HP/Documents/R/win-library/4.0/pomp/include/pomp.h>
#include <R_ext/Rdynload.h>

 


/* C snippet: 'rinit' */
#define Beta1		(__p[__parindex[0]])
#define Beta2		(__p[__parindex[1]])
#define mu_S1S2		(__p[__parindex[2]])
#define mu_IR		(__p[__parindex[3]])
#define mu_RS1		(__p[__parindex[4]])
#define eta		(__p[__parindex[5]])
#define rho		(__p[__parindex[6]])
#define N		(__p[__parindex[7]])
#define S1		(__x[__stateindex[0]])
#define S2		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define H		(__x[__stateindex[4]])

void __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)
{
 S1 = nearbyint(eta*N);
                      S2 = 0;
                      I = 1;
                      R = 0; 
                      H = 0; 
}

#undef Beta1
#undef Beta2
#undef mu_S1S2
#undef mu_IR
#undef mu_RS1
#undef eta
#undef rho
#undef N
#undef S1
#undef S2
#undef I
#undef R
#undef H

/* C snippet: 'step.fn' */
#define Beta1		(__p[__parindex[0]])
#define Beta2		(__p[__parindex[1]])
#define mu_S1S2		(__p[__parindex[2]])
#define mu_IR		(__p[__parindex[3]])
#define mu_RS1		(__p[__parindex[4]])
#define eta		(__p[__parindex[5]])
#define rho		(__p[__parindex[6]])
#define N		(__p[__parindex[7]])
#define S1		(__x[__stateindex[0]])
#define S2		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define H		(__x[__stateindex[4]])

void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t, double dt)
{
 double dN_S1S2 = rbinom(S1, 1-exp(-mu_S1S2*dt));
                      double dN_S1I = rbinom(S1 - dN_S1S2, 1-exp(-Beta1*I/N*dt));
                      double dN_S2I = rbinom(S2, 1-exp(-Beta2*I/N*dt));
                      double dN_IR = rbinom(I, 1-exp(-mu_IR*dt));
                      double dN_RS1 = rbinom(R, 1-exp(-mu_RS1*dt));
                      S1 += dN_RS1 - (dN_S1I + dN_S1S2);
                      S2 += dN_S1S2 - dN_S2I;
                      I += dN_S1I + dN_S2I - dN_IR;
                      R += dN_IR - dN_RS1;
                      H += dN_S1I + dN_S2I; 
}

#undef Beta1
#undef Beta2
#undef mu_S1S2
#undef mu_IR
#undef mu_RS1
#undef eta
#undef rho
#undef N
#undef S1
#undef S2
#undef I
#undef R
#undef H

/* C snippet: 'rmeasure' */
#define Beta1		(__p[__parindex[0]])
#define Beta2		(__p[__parindex[1]])
#define mu_S1S2		(__p[__parindex[2]])
#define mu_IR		(__p[__parindex[3]])
#define mu_RS1		(__p[__parindex[4]])
#define eta		(__p[__parindex[5]])
#define rho		(__p[__parindex[6]])
#define N		(__p[__parindex[7]])
#define S1		(__x[__stateindex[0]])
#define S2		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define H		(__x[__stateindex[4]])
#define count		(__y[__obsindex[0]])

void __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 count = rbinom(H,rho); 
}

#undef Beta1
#undef Beta2
#undef mu_S1S2
#undef mu_IR
#undef mu_RS1
#undef eta
#undef rho
#undef N
#undef S1
#undef S2
#undef I
#undef R
#undef H
#undef count

/* C snippet: 'dmeasure' */
#define Beta1		(__p[__parindex[0]])
#define Beta2		(__p[__parindex[1]])
#define mu_S1S2		(__p[__parindex[2]])
#define mu_IR		(__p[__parindex[3]])
#define mu_RS1		(__p[__parindex[4]])
#define eta		(__p[__parindex[5]])
#define rho		(__p[__parindex[6]])
#define N		(__p[__parindex[7]])
#define S1		(__x[__stateindex[0]])
#define S2		(__x[__stateindex[1]])
#define I		(__x[__stateindex[2]])
#define R		(__x[__stateindex[3]])
#define H		(__x[__stateindex[4]])
#define count		(__y[__obsindex[0]])
#define lik		(__lik[0])

void __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars, double t)
{
 lik = dbinom(count,H,rho,give_log);
                  
}

#undef Beta1
#undef Beta2
#undef mu_S1S2
#undef mu_IR
#undef mu_RS1
#undef eta
#undef rho
#undef N
#undef S1
#undef S2
#undef I
#undef R
#undef H
#undef count
#undef lik

/* C snippet: 'toEst' */
#define Beta1		(__p[__parindex[0]])
#define Beta2		(__p[__parindex[1]])
#define mu_S1S2		(__p[__parindex[2]])
#define mu_IR		(__p[__parindex[3]])
#define mu_RS1		(__p[__parindex[4]])
#define eta		(__p[__parindex[5]])
#define rho		(__p[__parindex[6]])
#define N		(__p[__parindex[7]])
#define T_Beta1		(__pt[__parindex[0]])
#define T_Beta2		(__pt[__parindex[1]])
#define T_mu_S1S2		(__pt[__parindex[2]])
#define T_mu_IR		(__pt[__parindex[3]])
#define T_mu_RS1		(__pt[__parindex[4]])
#define T_eta		(__pt[__parindex[5]])
#define T_rho		(__pt[__parindex[6]])
#define T_N		(__pt[__parindex[7]])

void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex)
{
 	T_N = log(N);
	T_Beta1 = log(Beta1);
	T_Beta2 = log(Beta2);
	T_mu_S1S2 = log(mu_S1S2);
	T_mu_IR = log(mu_IR);
	T_mu_RS1 = log(mu_RS1);
	T_rho = logit(rho);
	T_eta = logit(eta); 
}

#undef Beta1
#undef Beta2
#undef mu_S1S2
#undef mu_IR
#undef mu_RS1
#undef eta
#undef rho
#undef N
#undef T_Beta1
#undef T_Beta2
#undef T_mu_S1S2
#undef T_mu_IR
#undef T_mu_RS1
#undef T_eta
#undef T_rho
#undef T_N

/* C snippet: 'fromEst' */
#define Beta1		(__p[__parindex[0]])
#define Beta2		(__p[__parindex[1]])
#define mu_S1S2		(__p[__parindex[2]])
#define mu_IR		(__p[__parindex[3]])
#define mu_RS1		(__p[__parindex[4]])
#define eta		(__p[__parindex[5]])
#define rho		(__p[__parindex[6]])
#define N		(__p[__parindex[7]])
#define T_Beta1		(__pt[__parindex[0]])
#define T_Beta2		(__pt[__parindex[1]])
#define T_mu_S1S2		(__pt[__parindex[2]])
#define T_mu_IR		(__pt[__parindex[3]])
#define T_mu_RS1		(__pt[__parindex[4]])
#define T_eta		(__pt[__parindex[5]])
#define T_rho		(__pt[__parindex[6]])
#define T_N		(__pt[__parindex[7]])

void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex)
{
 	N = exp(T_N);
	Beta1 = exp(T_Beta1);
	Beta2 = exp(T_Beta2);
	mu_S1S2 = exp(T_mu_S1S2);
	mu_IR = exp(T_mu_IR);
	mu_RS1 = exp(T_mu_RS1);
	rho = expit(T_rho);
	eta = expit(T_eta); 
}

#undef Beta1
#undef Beta2
#undef mu_S1S2
#undef mu_IR
#undef mu_RS1
#undef eta
#undef rho
#undef N
#undef T_Beta1
#undef T_Beta2
#undef T_mu_S1S2
#undef T_mu_IR
#undef T_mu_RS1
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
