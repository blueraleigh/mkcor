#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

SEXP C_ephylo(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_ephylo_preorder(SEXP,SEXP,SEXP,SEXP);
SEXP C_ephylo_preorder_tips(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_ephylo_postorder(SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_ephylo_read_newick(SEXP,SEXP);

SEXP C_mkcor_em(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,
    SEXP,SEXP,SEXP);
SEXP C_mkcor_loglik(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
SEXP C_mkcor_simulate(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(C_ephylo, 5),
    CALLDEF(C_ephylo_preorder, 4),
    CALLDEF(C_ephylo_preorder_tips, 5),
    CALLDEF(C_ephylo_postorder, 5),
    CALLDEF(C_ephylo_read_newick, 2),
    CALLDEF(C_mkcor_em, 15),
    CALLDEF(C_mkcor_loglik, 11),
    CALLDEF(C_mkcor_simulate, 9),
    {NULL, NULL, 0}
};

void attribute_visible R_init_mkcor(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
