#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <assert.h>

#ifndef __GNUC__
#define  __attribute__(x) /* empty */
#elif __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 5) || __STRICT_ANSI__
#define  __attribute__(x)
#endif

#define SQRT_DBL_EPSILON 1.490116119384765696e-8

// test for a == b
// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
int
fequals(double a, double b)
{
    if (!R_FINITE(a) || !R_FINITE(b)) return 0;
    double diff = fabs(a - b);
    a = fabs(a);
    b = fabs(b);
    double M = (b > a) ? b : a;
    // use relative tolerance, modelled after R's all.equal.numeric
    if (M > SQRT_DBL_EPSILON)
        return (diff <= M * SQRT_DBL_EPSILON) ? 1 : 0;
    // otherwise, use absolute tolerance
    return (diff <= SQRT_DBL_EPSILON) ? 1 : 0;
}


enum {
    P00,
    P10,
    P01,
    P11
};

enum {
    P00_P00,
    P00_P10,
    P00_P01,
    P00_P11,
    P10_P10,
    P10_P01,
    P10_P11,
    P01_P01,
    P01_P11,
    P11_P11
};

static int CONVOLVE_MAP[4][4] = {
    {P00_P00, P00_P10, P00_P01, P00_P11},
    {P00_P10, P10_P10, P10_P01, P10_P11},
    {P00_P01, P10_P01, P01_P01, P01_P11},
    {P00_P11, P10_P11, P01_P11, P11_P11}
};

typedef struct ctmc {
    /* number of states in the model */
    int num_states;
    int num_states_x;
    int num_states_y;
    /* overall rates (rx, ry, rxy) */
    double rates[3];
    /* transition rates */
    double trates[3];
    /* transition probabilities */
    double pij[4];
    double log_pij[4];
    /* space for computing local expectations */
    double coef[10][16];
    double convolve[10];
    double I[16];
    /* non-zero eigenvalues */
    double values[3];
    /* index the transition type */
    int *tindex;
    /* space for log summation */
    double *work;
} ctmc_t;


static void
EXPAND(double p[4], double q[4], double A[16])
{
    A[0] = p[0]*q[0];
    A[1] = p[1]*q[1];
    A[2] = p[2]*q[2];
    A[3] = p[3]*q[3];
    A[4] = p[1]*q[0];
    A[5] = p[0]*q[1];
    A[6] = p[2]*q[0];
    A[7] = p[0]*q[2];
    A[8] = p[3]*q[0];
    A[9] = p[0]*q[3];
    A[10] = p[1]*q[2];
    A[11] = p[2]*q[1];
    A[12] = p[1]*q[3];
    A[13] = p[3]*q[1];
    A[14] = p[2]*q[3];
    A[15] = p[3]*q[2];
}


static void
ctmc_expand_coefficients(ctmc_t *ctmc)
{
    double kx = ctmc->num_states_x;
    double ky = ctmc->num_states_y;
    double p00[4] = {1, (kx-1), (ky-1), (kx-1)*(ky-1)};
    double p10[4] = {1, -1, (ky-1), -(ky-1)};
    double p01[4] = {1, (kx-1), -1, -(kx-1)};
    double p11[4] = {1, -1, -1, 1};
    EXPAND(p00, p00, ctmc->coef[P00_P00]);
    EXPAND(p00, p10, ctmc->coef[P00_P10]);
    EXPAND(p00, p01, ctmc->coef[P00_P01]);
    EXPAND(p00, p11, ctmc->coef[P00_P11]);
    EXPAND(p10, p10, ctmc->coef[P10_P10]);
    EXPAND(p10, p01, ctmc->coef[P10_P01]);
    EXPAND(p10, p11, ctmc->coef[P10_P11]);
    EXPAND(p01, p01, ctmc->coef[P01_P01]);
    EXPAND(p01, p11, ctmc->coef[P01_P11]);
    EXPAND(p11, p11, ctmc->coef[P11_P11]);
}


static void
ctmc_set_transition_type(ctmc_t *ctmc, int a, int b)
{
    int num_states_y = ctmc->num_states_y;
    int k = a % num_states_y;
    int i = (a - k) / num_states_y;
    int l = b % num_states_y;
    int j = (b - l) / num_states_y;
    int i1 = a + b*ctmc->num_states;
    int i2 = b + a*ctmc->num_states;
    if (i == j && k == l)
    {
        ctmc->tindex[i1] = ctmc->tindex[i2] = P00;
    }
    else if (i != j && k == l)
    {
        ctmc->tindex[i1] = ctmc->tindex[i2] = P10;
    }
    else if (i == j && k != l)
    {
        ctmc->tindex[i1] = ctmc->tindex[i2] = P01;
    }
    else
    {
        ctmc->tindex[i1] = ctmc->tindex[i2] = P11;
    }
}


static int
ctmc_get_transition_type(ctmc_t *ctmc, int a, int b)
{
    return ctmc->tindex[a+b*ctmc->num_states];
}


static void
ctmc_init(ctmc_t *ctmc, int num_states_x, int num_states_y)
{
    int i;
    int j;
    int num_states = num_states_x * num_states_y;
    ctmc->num_states = num_states;
    ctmc->num_states_x = num_states_x;
    ctmc->num_states_y = num_states_y;
    ctmc->tindex = (int *)R_alloc(num_states * num_states, sizeof(int));
    ctmc->work = (double *)R_alloc(num_states, sizeof(double));
    for (i = 0; i < num_states; ++i)
    {
        for (j = i; j < num_states; ++j)
            ctmc_set_transition_type(ctmc, i, j);
    }
    ctmc_expand_coefficients(ctmc);
}


/* Eigenvectors of the transition rate matrix (currently unused) */
__attribute__((unused)) static void
ctmc_eigenvectors(ctmc_t *ctmc, double *vectors)
{
    int i;
    int j;
    int k;
    int num_states = ctmc->num_states;
    int num_states_x = ctmc->num_states_x;
    int num_states_y = ctmc->num_states_y;
    double *vx = vectors + num_states;
    double *vy = vx + num_states*(num_states_x-1);
    double *vxy = vy + num_states*(num_states_y-1);
    Memzero(vectors, num_states*num_states);
    for (i = 0; i < num_states; ++i)
        vectors[i] = 1;
    for (j = 0; j < num_states_x-1; ++j)
    {
        for (i = 0; i < num_states_y; ++i)
            vx[i+j*num_states] = 1;
        for (i = (j+1)*num_states_y; i < 2*(j+1)*num_states_y; ++i)
            vx[i+j*num_states] = -1;
    }
    for (j = 0; j < num_states_y-1; ++j)
    {
        for (i = 0; i < num_states; i += num_states_y)
            vy[i+j*num_states] = 1;
        for (i = (j+1); i < num_states; i += num_states_y)
            vy[i+j*num_states] = -1;
    }
    for (j = 0; j < (num_states_x-1)*(num_states_y-1); ++j)
    {
        vxy[0+j*num_states] = 1;
        vxy[(j+1)*(num_states_y+1)+j*num_states] = 1;
        vxy[(j+1)+j*num_states] = -1;
        vxy[(j+1)*(num_states_y)+j*num_states] = -1;
    }
    // orthogonalize the eigenvectors with modified gram-schmidt
    int one = 1;
    double scal;
    double *a;
    double *b;
    for (i = 1; i < num_states; ++i)
    {
        b = vectors + i*num_states;
        for (j = 0; j < i; ++j)
        {
            a = vectors + j*num_states;
            scal = F77_CALL(ddot)(&num_states, a, &one, b, &one) /
                F77_CALL(ddot)(&num_states, a, &one, a, &one);
            for (k = 0; k < num_states; ++k)
                b[k] -= scal*a[k];
        }
    }
    // normalize the eigenvectors
    for (i = 0; i < num_states; ++i)
    {
        a = vectors + i*num_states;
        scal = 1 / F77_CALL(dnrm2)(&num_states, a, &one);
        F77_CALL(dscal)(&num_states, &scal, a, &one);
    }
}


static void
ctmc_set_rates(ctmc_t *ctmc, double *rates)
{
    double kx = ctmc->num_states_x;
    double ky = ctmc->num_states_y;
    ctmc->rates[0] = rates[0];
    ctmc->rates[1] = rates[1];
    ctmc->rates[2] = rates[2];
    ctmc->trates[0] = rates[0] / (kx - 1);
    ctmc->trates[1] = rates[1] / (ky - 1);
    ctmc->trates[2] = rates[2] / ((kx - 1)*(ky - 1));
    ctmc->values[0] = -kx*(rates[0] + rates[2])/(kx-1);
    ctmc->values[1] = -ky*(rates[1] + rates[2])/(ky-1);
    ctmc->values[2] = -(kx*rates[0]/(kx-1) + ky*rates[1]/(ky-1) + 
        ((kx-1)*(ky-1)-1)*rates[2]/((kx-1)*(ky-1)));
}


static void
ctmc_set_log_probabilities(ctmc_t *ctmc, double time)
{
    double k = ctmc->num_states;
    double kx = ctmc->num_states_x;
    double ky = ctmc->num_states_y;
    double dx = ctmc->values[0]*time;
    double dy = ctmc->values[1]*time;
    double dxy = ctmc->values[2]*time;
    double log_kx = log(kx-1);
    double log_ky = log(ky-1);
    double log_kxy = log_kx + log_ky;
    double log_k = log(k);
    ctmc->log_pij[P00] = Rf_logspace_add(
        Rf_logspace_add(0, log_kx + dx), 
        Rf_logspace_add(log_ky + dy, log_kxy + dxy)) - log_k;
    ctmc->log_pij[P10] = Rf_logspace_sub(
        Rf_logspace_add(0, log_ky + dy), 
        Rf_logspace_add(dx, log_ky + dxy)) - log_k;
    ctmc->log_pij[P01] = Rf_logspace_sub(
        Rf_logspace_add(0, log_kx + dx), 
        Rf_logspace_add(dy, log_kx + dxy)) - log_k;
    ctmc->log_pij[P11] = Rf_logspace_sub(
        Rf_logspace_add(0, dxy), 
        Rf_logspace_add(dx, dy)) - log_k;
    // fix possible NaN's resulting from -Inf - (-Inf)
    if (ISNAN(ctmc->log_pij[P00])) ctmc->log_pij[P00] = R_NegInf;
    if (ISNAN(ctmc->log_pij[P10])) ctmc->log_pij[P10] = R_NegInf;
    if (ISNAN(ctmc->log_pij[P01])) ctmc->log_pij[P01] = R_NegInf;
    if (ISNAN(ctmc->log_pij[P11])) ctmc->log_pij[P11] = R_NegInf;
}


static void
ctmc_set_probabilities(ctmc_t *ctmc, double time)
{
    double k = ctmc->num_states;
    double kx = ctmc->num_states_x;
    double ky = ctmc->num_states_y;
    double Dx = exp(ctmc->values[0]*time);
    double Dy = exp(ctmc->values[1]*time);
    double Dxy = exp(ctmc->values[2]*time);
    ctmc->pij[P00] = (1 + (kx-1)*Dx + (ky-1)*Dy + (kx-1)*(ky-1)*Dxy) / k;
    ctmc->pij[P10] = (1 - Dx + (ky-1)*Dy - (ky-1)*Dxy) / k;
    ctmc->pij[P01] = (1 + (kx-1)*Dx - Dy - (kx-1)*Dxy) / k;
    ctmc->pij[P11] = (1 - Dx - Dy + Dxy) / k;
    // fix negative 0
    ctmc->pij[P00] = ctmc->pij[P00] <= 0.0 ? 0.0 : ctmc->pij[P00];
    ctmc->pij[P10] = ctmc->pij[P10] <= 0.0 ? 0.0 : ctmc->pij[P10];
    ctmc->pij[P01] = ctmc->pij[P01] <= 0.0 ? 0.0 : ctmc->pij[P01];
    ctmc->pij[P11] = ctmc->pij[P11] <= 0.0 ? 0.0 : ctmc->pij[P11];
}


static double
ctmc_get_probability(ctmc_t *ctmc, int a, int b)
{
    return ctmc->pij[ctmc_get_transition_type(ctmc, a, b)];
}


static double
ctmc_get_log_probability(ctmc_t *ctmc, int a, int b)
{
    return ctmc->log_pij[ctmc_get_transition_type(ctmc, a, b)];
}


static double
ctmc_get_rate(ctmc_t *ctmc, int a, int b)
{
    if (a == b) return 0;
    return ctmc->trates[ctmc_get_transition_type(ctmc, a, b) - 1];
}


static double
CONVOLVE(double A[16], double B[16], double time)
{
    double c = time + 
        A[1]*B[1] + 
        A[2]*B[2] +
        A[3]*B[3] +
        A[4]*B[4] + 
        A[5]*B[5] + 
        A[6]*B[6] +
        A[7]*B[7] +
        A[8]*B[8] +
        A[9]*B[9] +
        A[10]*B[10] +
        A[11]*B[11] +
        A[12]*B[12] +
        A[13]*B[13] +
        A[14]*B[14] +
        A[15]*B[15];
    // fix negative 0
    return (c <= 0) ? 0.0 : c;
}


static void
ctmc_convolve(ctmc_t *ctmc, double time)
{
    double dx = ctmc->values[0];
    double dy = ctmc->values[1];
    double dxy = ctmc->values[2];
    double Dx = exp(dx*time);
    double Dy = exp(dy*time);
    double Dxy = exp(dxy*time);
    double k = ctmc->num_states * ctmc->num_states;
    ctmc->I[0] = time * 1;
    ctmc->I[1] = time * Dx;
    ctmc->I[2] = time * Dy;
    ctmc->I[3] = time * Dxy;
    if (fequals(dx, 0))
    {
        ctmc->I[4] = time * Dx;
        ctmc->I[5] = time * 1;
    }
    else
    {
        ctmc->I[4] = (Dx-1)/dx;
        ctmc->I[5] = -(1-Dx)/dx;
    }
    if (fequals(dy, 0))
    {
        ctmc->I[6] = time * Dy;
        ctmc->I[7] = time * 1;
    }
    else
    {
        ctmc->I[6] = (Dy-1)/dy;
        ctmc->I[7] = -(1-Dy)/dy;
    }
    if (fequals(dxy, 0))
    {
        ctmc->I[8] = time * Dxy;
        ctmc->I[9] = time * 1;
    }
    else
    {
        ctmc->I[8] = (Dxy-1)/dxy;
        ctmc->I[9] = -(1-Dxy)/dxy;
    }
    if (fequals(dx, dy))
    {
        ctmc->I[10] = time * Dx;
        ctmc->I[11] = time * Dy;
    }
    else
    {
        ctmc->I[10] = (Dx-Dy)/(dx-dy);
        ctmc->I[11] = (Dy-Dx)/(dy-dx);
    }
    if (fequals(dx, dxy))
    {
        ctmc->I[12] = time * Dx;
        ctmc->I[13] = time * Dxy;
    }
    else
    {
        ctmc->I[12] = (Dx-Dxy)/(dx-dxy);
        ctmc->I[13] = (Dxy-Dx)/(dxy-dx);
    }
    if (fequals(dy, dxy))
    {
        ctmc->I[14] = time * Dy;
        ctmc->I[15] = time * Dxy;
    }
    else
    {
        ctmc->I[14] = (Dy-Dxy)/(dy-dxy);
        ctmc->I[15] = (Dxy-Dy)/(dxy-dy);
    }
    /* Integrate P00(t)*P00(T - t)dt from 0 to T */
    ctmc->convolve[P00_P00] = CONVOLVE(ctmc->coef[P00_P00], ctmc->I, time) / k;
    /* Integrate P00(t)*P10(T - t)dt from 0 to T */
    ctmc->convolve[P00_P10] = CONVOLVE(ctmc->coef[P00_P10], ctmc->I, time) / k;
    /* Integrate P00(t)*P01(T - t)dt from 0 to T */
    ctmc->convolve[P00_P01] = CONVOLVE(ctmc->coef[P00_P01], ctmc->I, time) / k;
    /* Integrate P00(t)*P11(T - t)dt from 0 to T */
    ctmc->convolve[P00_P11] = CONVOLVE(ctmc->coef[P00_P11], ctmc->I, time) / k;
    /* Integrate P10(t)*P10(T - t)dt from 0 to T */
    ctmc->convolve[P10_P10] = CONVOLVE(ctmc->coef[P10_P10], ctmc->I, time) / k;
    /* Integrate P10(t)*P01(T - t)dt from 0 to T */
    ctmc->convolve[P10_P01] = CONVOLVE(ctmc->coef[P10_P01], ctmc->I, time) / k;
    /* Integrate P10(t)*P11(T - t)dt from 0 to T */
    ctmc->convolve[P10_P11] = CONVOLVE(ctmc->coef[P10_P11], ctmc->I, time) / k;
    /* Integrate P01(t)*P01(T - t)dt from 0 to T */
    ctmc->convolve[P01_P01] = CONVOLVE(ctmc->coef[P01_P01], ctmc->I, time) / k;
    /* Integrate P01(t)*P11(T - t)dt from 0 to T */
    ctmc->convolve[P01_P11] = CONVOLVE(ctmc->coef[P01_P11], ctmc->I, time) / k;
    /* Integrate P11(t)*P11(T - t)dt from 0 to T */
    ctmc->convolve[P11_P11] = CONVOLVE(ctmc->coef[P11_P11], ctmc->I, time) / k;
}

/* Compute the expected number of i -> j transitions on a branch of
** length time given that it started in state a and ended in state b.
**
** See appendix in:
** Hobolth, A., & Jensen, J. L. (2005). Statistical inference in
** evolutionary models of DNA sequences via the EM algorithm.
** Statistical applications in genetics and molecular biology, 4(1).
*/
static double
ctmc_get_local_expectation(ctmc_t *ctmc, int a, int i, int j, int b)
{
    if (i == j) return 0;
    double qij = ctmc_get_rate(ctmc, i, j);
    if (fequals(qij, 0)) return 0;
    //double pab = ctmc_get_probability(ctmc, a, b);
    //if (fequals(pab, 0)) return 0;
    double log_pab = ctmc_get_log_probability(ctmc, a, b);
    if (!R_FINITE(log_pab)) return 0;
    int tai = ctmc_get_transition_type(ctmc, a, i);
    int tjb = ctmc_get_transition_type(ctmc, j, b);
    //return qij * ctmc->convolve[CONVOLVE_MAP[tai][tjb]] / pab;
    return exp(log(qij) + log(ctmc->convolve[CONVOLVE_MAP[tai][tjb]]) - log_pab);
}


/*
Schematic of how inside and outside likelihoods relate to the tree
structure. In the example, these likelihoods apply to the internal
node marked by `.`

D = `inside likelihood crown`  : likelihood of data in subtree to A and B
S = `inside likelihood stem`   : likelihood of data in subtree to A and B
F = `outside likelihood crown` : likelihood of data in subtree to C and D
U = `outside likelihood stem`  : likelihood of data in subtree to C and D

   A        
    \       B
     \     / 
      \   /     
       \ /     C
D,U --> .     /
         \   /
          \ /       
   S,F --> +       D
            \     / 
             \   /
              \ /
               +
*/


static void
marginal_ancestral_state_prob(
    int v,
    double logL,
    double *inside_loglik_crown,
    double *outside_loglik_stem,
    double *state_prob,
    ctmc_t *ctmc)
{
    int i;
    int num_states = ctmc->num_states;
    int vn = v*num_states;

    double *D = inside_loglik_crown + vn;
    double *U = outside_loglik_stem + vn;
    double *Z = state_prob + vn;

    for (i = 0; i < num_states; ++i)
        Z[i] = exp(D[i] + U[i] - logL);
}


static void
inside_loglikelihood_stem(
    int u,
    double brlen,
    double *inside_loglik_stem,
    double *inside_loglik_crown,
    ctmc_t *ctmc)
{
    int i;
    int j;
    int num_states = ctmc->num_states;
    int offset = u*num_states;
    double *tmp = ctmc->work;
    double *D = inside_loglik_crown + offset;
    double *S = inside_loglik_stem + offset;

    //ctmc_set_probabilities(ctmc, brlen);
    ctmc_set_log_probabilities(ctmc, brlen);

    for (i = 0; i < num_states; ++i)
    {
        for (j = 0; j < num_states; ++j)
        {
            //tmp[j] = log(ctmc_get_probability(ctmc, i, j)) + D[j];
            tmp[j] = ctmc_get_log_probability(ctmc, i, j) + D[j];
        }
        S[i] = Rf_logspace_sum(tmp, num_states);
        if (ISNAN(S[i])) S[i] = R_NegInf;
    }
}


static void
inside_loglikelihood_crown(
    int u,
    int *left_child,
    int *right_sib,
    double *inside_loglik_stem,
    double *inside_loglik_crown,
    ctmc_t *ctmc)
{
    int k;
    int v;
    int num_states = ctmc->num_states;
    double *D;
    double *S;

    D = inside_loglik_crown + u*num_states;

    Memzero(D, num_states);

    for (v = left_child[u]-1; v != -1; v = right_sib[v]-1)
    {
        S = inside_loglik_stem + v*num_states;
        for (k = 0; k < num_states; ++k)
            D[k] += S[k];
    }
}


static void
outside_loglikelihood_crown(
    int v,
    int *left_child,
    int *right_sib,
    int *parent,
    double *inside_loglik_stem,
    double *outside_loglik_crown,
    double *outside_loglik_stem,
    ctmc_t *ctmc)
{
    int k;
    int s;
    int num_states = ctmc->num_states;
    int u = parent[v] - 1;
    int un = u*num_states;
    double p = -log(num_states);

    double *S;
    double *F = outside_loglik_crown;

    if (u == -1)
    {
        F = outside_loglik_stem + v*num_states;
        for (k = 0; k < num_states; ++k)
            F[k] = p;
        return;
    }
    for (k = 0; k < num_states; ++k)
        F[k] = outside_loglik_stem[k + un];
    for (s = left_child[u]-1; s != -1; s = right_sib[s]-1)
    {
        if (s == v) continue;
        S = inside_loglik_stem + s*num_states;
        for (k = 0; k < num_states; ++k)
            F[k] += S[k];
    }
}


static void
outside_loglikelihood_stem(
    int v,
    double brlen,
    double logL,
    double *inside_loglik_crown,
    double *outside_loglik_crown,
    double *outside_loglik_stem,
    double *state_prob,
    double *branch_counts,
    double *expected_counts,
    ctmc_t *ctmc)
{
    int i;
    int j;
    int k;
    int l;
    int tindex;
    int num_states = ctmc->num_states;
    int num_states2 = num_states * num_states;
    int vn = v*num_states;

    double e1;
    double e2;
    double P;
    double logP;
    double *tmp = ctmc->work;

    double *F = outside_loglik_crown;
    double *D = inside_loglik_crown + vn;
    double *U = outside_loglik_stem + vn;
    double *Z = state_prob + vn;
    double *E = branch_counts + num_states2*v;
    double *N = expected_counts;

    //ctmc_set_probabilities(ctmc, brlen);
    ctmc_set_log_probabilities(ctmc, brlen);
    ctmc_convolve(ctmc, brlen);

    Memzero(E, num_states2);

    for (i = 0; i < num_states; ++i)
    {
        Z[i] = 0;
        for (j = 0; j < num_states; ++j)
        {
            //tmp[j] = F[j] + log(ctmc_get_probability(ctmc, j, i));
            tmp[j] = F[j] + ctmc_get_log_probability(ctmc, j, i);
            logP = tmp[j] + D[i];
            // posterior log probability of j -> i transition on branch
            P = exp(logP - logL);
            assert(!ISNAN(P));
            Z[i] += P;
            for (k = 0; k < num_states; ++k)
            {
                for (l = k+1; l < num_states; ++l)
                {
                    tindex = ctmc_get_transition_type(ctmc, k, l) - 1;
                    e1 = P * ctmc_get_local_expectation(ctmc, j, k, l, i);
                    e2 = P * ctmc_get_local_expectation(ctmc, j, l, k, i);
                    assert(!ISNAN(e1));
                    assert(!ISNAN(e2));
                    N[tindex] += e1 + e2;
                    E[k+l*num_states] += e1;
                    E[l+k*num_states] += e2;
                }
            }
        }
        U[i] = Rf_logspace_sum(tmp, num_states);
        if (ISNAN(U[i])) U[i] = R_NegInf;
    }
}


// log likelihood of the data in the subtree rooted at node u
static double
loglikelihood(int u, double *inside_loglik_crown, ctmc_t *ctmc)
{
    int num_states = ctmc->num_states;
    double *D = inside_loglik_crown + u*num_states;
    double logL = Rf_logspace_sum(D, num_states) - log(num_states);
    if (ISNAN(logL)) Rf_error("non-finite likelihood");
    return logL;
}


static double
estep(
    int num_tips,
    int num_nodes,
    int *postorder,
    int *preorder,
    int *parent,
    int *left_child,
    int *right_sib,
    double *brlen,
    double *inside_loglik_stem,
    double *inside_loglik_crown,
    double *outside_loglik_stem,
    double *outside_loglik_crown,
    double *state_prob,
    double *branch_counts,
    double *expected_counts,
    ctmc_t *ctmc)
{
    int i;
    int u;
    int v;
    double logL;
    Memzero(expected_counts, 3);
    for (i = 0; i < num_nodes; ++i)
    {
        u = postorder[i] - 1;
        if (u >= num_tips)
        {
            inside_loglikelihood_crown(u, left_child, right_sib,
                inside_loglik_stem, inside_loglik_crown, ctmc);
        }
        if (u != num_tips)
        {
            inside_loglikelihood_stem(u, brlen[u], inside_loglik_stem, 
                inside_loglik_crown, ctmc);
        }
    }

    logL = loglikelihood(num_tips, inside_loglik_crown, ctmc);

    for (i = 0; i < num_nodes; ++i)
    {
        v = preorder[i] - 1;

        outside_loglikelihood_crown(v, left_child, right_sib, 
            parent, inside_loglik_stem, outside_loglik_crown, 
            outside_loglik_stem, ctmc);

        if (v == num_tips) 
        {
            // get the state probabilities for the root
            marginal_ancestral_state_prob(v, logL, inside_loglik_crown, 
                outside_loglik_stem, state_prob, ctmc);
            continue;
        }

        outside_loglikelihood_stem(v, brlen[v], logL, inside_loglik_crown,
            outside_loglik_crown, outside_loglik_stem, state_prob,
            branch_counts, expected_counts, ctmc);
    }
    return logL;
}


static void
mstep(
    double treelen,
    double *par,
    double *expected_counts,
    ctmc_t *ctmc)
{
    par[0] = expected_counts[0] / treelen;
    par[1] = expected_counts[1] / treelen;
    par[2] = expected_counts[2] / treelen;
    ctmc_set_rates(ctmc, par);
}


SEXP
C_mkcor_em(
    SEXP r_num_states,
    SEXP r_num_states_x,
    SEXP r_num_states_y,
    SEXP r_par,
    SEXP r_tol,
    SEXP r_maxiter,
    SEXP r_node_conditional_loglik,
    SEXP r_num_tips,
    SEXP r_num_nodes,
    SEXP r_postorder,
    SEXP r_preorder,
    SEXP r_parent,
    SEXP r_left_child,
    SEXP r_right_sib,
    SEXP r_brlen)
{
    int num_tips = *INTEGER(r_num_tips);
    int num_nodes = *INTEGER(r_num_nodes);
    int *postorder = INTEGER(r_postorder);
    int *preorder = INTEGER(r_preorder);
    int *parent = INTEGER(r_parent);
    int *left_child = INTEGER(r_left_child);
    int *right_sib = INTEGER(r_right_sib);
    double *brlen = REAL(r_brlen);

    int num_states = *INTEGER(r_num_states);
    int num_states_x = *INTEGER(r_num_states_x);
    int num_states_y = *INTEGER(r_num_states_y);

    double *inside_loglik_crown = REAL(r_node_conditional_loglik);
    double *inside_loglik_stem = (double *) R_alloc(
        num_states*num_nodes, sizeof(*inside_loglik_stem));
    double *outside_loglik_stem = (double *) R_alloc(
        num_states*num_nodes, sizeof(*outside_loglik_stem));
    double *outside_loglik_crown = (double *) R_alloc(
        num_states, sizeof(*outside_loglik_crown));

    double par[3] = {0};
    Memcpy(par, REAL(r_par), Rf_length(r_par));

    ctmc_t ctmc;

    ctmc_init(&ctmc, num_states_x, num_states_y);
    ctmc_set_rates(&ctmc, par);

    int it;
    int maxit = *INTEGER(r_maxiter);

    SEXP ans1 = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP ans2 = PROTECT(Rf_allocVector(REALSXP, 3));
    SEXP ans3 = PROTECT(Rf_allocMatrix(REALSXP, num_states, num_nodes));
    SEXP dims = PROTECT(Rf_allocVector(INTSXP, 3));
    SET_INTEGER_ELT(dims, 0, num_states);
    SET_INTEGER_ELT(dims, 1, num_states);
    SET_INTEGER_ELT(dims, 2, num_nodes);
    SEXP ans4 = PROTECT(Rf_allocArray(REALSXP, dims));
    double *expected_counts = REAL(ans2);
    double *state_prob = REAL(ans3);
    double *branch_counts = REAL(ans4);

    Memzero(branch_counts, num_states*num_states*num_nodes);

    double tol = *REAL(r_tol);
    double delta;
    double logL;
    double logL0 = R_NegInf;

    double *mll = (double *) R_alloc(maxit, sizeof(*mll));

    double treelen = 0;
    for (it = 0; it < num_nodes; ++it)
        treelen += brlen[it];

    for (it = 0; it < maxit; ++it)
    {
        logL = estep(
            num_tips, 
            num_nodes,
            postorder,
            preorder,
            parent,
            left_child,
            right_sib,
            brlen,
            inside_loglik_stem,
            inside_loglik_crown,
            outside_loglik_stem, 
            outside_loglik_crown,
            state_prob,
            branch_counts,
            expected_counts,
            &ctmc
        );
        mll[it] = logL;
        delta = logL - logL0;
        if (delta < tol) break;
        mstep(
            treelen,
            par,
            expected_counts,
            &ctmc
        );
        logL0 = logL;
    }
    if (it == maxit) --it;
    SEXP ans5 = PROTECT(Rf_allocVector(REALSXP, it+1));
    Memcpy(REAL(ans1), par, 3);
    Memcpy(REAL(ans5), mll, it+1);
    UNPROTECT(6);
    return Rf_list5(ans1, ans2, ans3, ans4, ans5);
}


SEXP
C_mkcor_loglik(
    SEXP r_num_states,
    SEXP r_num_states_x,
    SEXP r_num_states_y,
    SEXP r_par,
    SEXP r_node_conditional_loglik,
    SEXP r_num_tips,
    SEXP r_num_nodes,
    SEXP r_postorder,
    SEXP r_left_child,
    SEXP r_right_sib,
    SEXP r_brlen)
{
    int num_tips = *INTEGER(r_num_tips);
    int num_nodes = *INTEGER(r_num_nodes);
    int *postorder = INTEGER(r_postorder);
    int *left_child = INTEGER(r_left_child);
    int *right_sib = INTEGER(r_right_sib);
    double *brlen = REAL(r_brlen);

    int num_states = *INTEGER(r_num_states);
    int num_states_x = *INTEGER(r_num_states_x);
    int num_states_y = *INTEGER(r_num_states_y);

    double *inside_loglik_crown = REAL(r_node_conditional_loglik);
    double *inside_loglik_stem = (double *) R_alloc(
        num_states*num_nodes, sizeof(*inside_loglik_stem));
    double *outside_loglik_stem = (double *) R_alloc(
        num_states*num_nodes, sizeof(*outside_loglik_stem));
    double *outside_loglik_crown = (double *) R_alloc(
        num_states, sizeof(*outside_loglik_crown));

    double par[3] = {0};
    Memcpy(par, REAL(r_par), Rf_length(r_par));

    ctmc_t ctmc;

    ctmc_init(&ctmc, num_states_x, num_states_y);
    ctmc_set_rates(&ctmc, par);

    int i;
    int u;
    double logL;

    for (i = 0; i < num_nodes; ++i)
    {
        u = postorder[i] - 1;
        if (u >= num_tips)
        {
            inside_loglikelihood_crown(u, left_child, right_sib,
                inside_loglik_stem, inside_loglik_crown, &ctmc);
        }
        if (u != num_tips)
        {
            inside_loglikelihood_stem(u, brlen[u], inside_loglik_stem, 
                inside_loglik_crown, &ctmc);
        }
    }

    logL = loglikelihood(num_tips, inside_loglik_crown, &ctmc);

    return Rf_ScalarReal(logL);
}


SEXP
C_mkcor_simulate(
    SEXP r_num_states,
    SEXP r_num_states_x,
    SEXP r_num_states_y,
    SEXP r_par,
    SEXP r_num_tips,
    SEXP r_num_nodes,
    SEXP r_preorder,
    SEXP r_parent,
    SEXP r_brlen)
{
    int a;
    int b;
    int x;
    int y;
    int i;
    int u;
    int v;
    int num_tips = *INTEGER(r_num_tips);
    int num_nodes = *INTEGER(r_num_nodes);
    int *preorder = INTEGER(r_preorder);
    int *parent = INTEGER(r_parent);
    double *brlen = REAL(r_brlen);

    int num_states = *INTEGER(r_num_states);
    int num_states_x = *INTEGER(r_num_states_x);
    int num_states_y = *INTEGER(r_num_states_y);
    
    double p;
    double par[3] = {0};
    Memcpy(par, REAL(r_par), Rf_length(r_par));

    ctmc_t ctmc;

    ctmc_init(&ctmc, num_states_x, num_states_y);
    ctmc_set_rates(&ctmc, par);

    int *node_state = (int *)R_alloc(num_nodes, sizeof(*node_state));
    SEXP r_node_state_x = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_node_state_y = PROTECT(Rf_allocVector(INTSXP, num_nodes));

    GetRNGstate();

    a = R_unif_index(num_states);
    y = a % num_states_y;
    x = (a - y) / num_states_y;
    node_state[num_tips] = a;
    SET_INTEGER_ELT(r_node_state_x, num_tips, x+1);
    SET_INTEGER_ELT(r_node_state_y, num_tips, y+1);
    for (i = 1; i < num_nodes; ++i)
    {
        v = preorder[i] - 1;
        u = parent[v] - 1;
        a = node_state[u];
        p = unif_rand();
        ctmc_set_probabilities(&ctmc, brlen[v]);
        for (b = 0; b < num_states; ++b)
        {
            p -= ctmc_get_probability(&ctmc, a, b);
            if (p <= 0) break;
        }
        y = b % num_states_y;
        x = (b - y) / num_states_y;
        node_state[v] = b;
        SET_INTEGER_ELT(r_node_state_x, v, x+1);
        SET_INTEGER_ELT(r_node_state_y, v, y+1);
    }

    PutRNGstate();

    UNPROTECT(2);
    return Rf_list2(r_node_state_x, r_node_state_y);
}
