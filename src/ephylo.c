#include <R.h>
#include <Rinternals.h>

SEXP
C_ephylo(
    SEXP r_num_nodes,
    SEXP r_num_edges,
    SEXP r_edge_parent,
    SEXP r_edge_child,
    SEXP r_edge_length)
{
    int num_nodes = *INTEGER(r_num_nodes);
    int num_edges = *INTEGER(r_num_edges);

    SEXP r_parent = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_left_child = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_right_child = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_left_sib = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_right_sib = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_brlen = PROTECT(Rf_allocVector(REALSXP, num_nodes));

    int *edge_parent = INTEGER(r_edge_parent);
    int *edge_child = INTEGER(r_edge_child);
    double *edge_length = REAL(r_edge_length);

    int *parent = INTEGER(r_parent);
    int *left_child = INTEGER(r_left_child);
    int *right_child = INTEGER(r_right_child);
    int *left_sib = INTEGER(r_left_sib);
    int *right_sib = INTEGER(r_right_sib);
    double *brlen = REAL(r_brlen);

    Memzero(parent, num_nodes);
    Memzero(left_child, num_nodes);
    Memzero(right_child, num_nodes);
    Memzero(left_sib, num_nodes);
    Memzero(right_sib, num_nodes);
    Memzero(brlen, num_nodes);

    int i;
    int u, u1;
    int v, v1;
    for (i = 0; i < num_edges; ++i)
    {
        u1 = edge_parent[i], u = u1 - 1;
        v1 = edge_child[i], v = v1 - 1;
        brlen[v] = edge_length[i];
        parent[v] = u1;
        if (left_child[u] == 0)
        {
            left_child[u] = v1;
        }
        else
        {
            if (right_child[u])
            {
                right_sib[right_child[u]-1] = v1;
                left_sib[v] = right_child[u];
            }
            else
            {
                right_sib[left_child[u]-1] = v1;
                left_sib[v] = left_child[u];
            }
            right_child[u] = v1;
        }
    }
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, r_parent);
    SET_VECTOR_ELT(ans, 1, r_left_child);
    SET_VECTOR_ELT(ans, 2, r_right_child);
    SET_VECTOR_ELT(ans, 3, r_left_sib);
    SET_VECTOR_ELT(ans, 4, r_right_sib);
    SET_VECTOR_ELT(ans, 5, r_brlen);
    UNPROTECT(7);
    return ans;
}

SEXP C_ephylo_preorder(
    SEXP r_num_nodes,
    SEXP r_root,
    SEXP r_right_child, 
    SEXP r_left_sib)
{
    int num_nodes = *INTEGER(r_num_nodes);
    int root = *INTEGER(r_root) - 1;
    int *right_child = INTEGER(r_right_child);
    int *left_sib = INTEGER(r_left_sib);

    int *stack = (int *) R_alloc(num_nodes, sizeof(*stack));
    
    int u;
    int v;
    int stack_top;

    int *preorder = (int *) R_alloc(num_nodes, sizeof(*preorder));

    stack_top = -1;
    for (u = (right_child[root]-1); u != -1; u = (left_sib[u]-1))
    {
        ++stack_top;
        stack[stack_top] = u;
    }

    int k = 1;
    while (stack_top >= 0)
    {
        u = stack[stack_top];
        preorder[k++] = u + 1;
        --stack_top;
        for (v = (right_child[u]-1); v != -1; v = (left_sib[v]-1))
        {
            ++stack_top;
            stack[stack_top] = v;
        }
    }
    *preorder = root + 1;
    
    SEXP r_preorder = PROTECT(Rf_allocVector(INTSXP, k));
    Memcpy(INTEGER(r_preorder), preorder, k);

    UNPROTECT(1);
    return r_preorder;
}

SEXP C_ephylo_postorder(
    SEXP r_num_nodes,
    SEXP r_root,
    SEXP r_parent,
    SEXP r_right_child,
    SEXP r_left_sib)
{
    int num_nodes = *INTEGER(r_num_nodes);
    int root = *INTEGER(r_root) - 1;
    int *parent = INTEGER(r_parent);
    int *right_child = INTEGER(r_right_child);
    int *left_sib = INTEGER(r_left_sib);

    int *stack = (int *) R_alloc(num_nodes, sizeof(*stack));
    
    int u;
    int v;
    int stack_top;

    int *postorder = (int *) R_alloc(num_nodes, sizeof(*postorder));

    stack_top = -1;
    for (u = (right_child[root]-1); u != -1; u = (left_sib[u]-1))
    {
        ++stack_top;
        stack[stack_top] = u;
    }

    int k = 0;
    int postorder_parent = -1;
    while (stack_top >= 0)
    {
        u = stack[stack_top];
        if ((u == postorder_parent) || (right_child[u] == 0))
        {
            --stack_top;
            postorder[k++] = u + 1;
            postorder_parent = parent[u] - 1;
        }
        else
        {
            for (v = (right_child[u]-1); v != -1; v = (left_sib[v]-1))
            {
                ++stack_top;
                stack[stack_top] = v;
            }
        }
    }
    postorder[k++] = postorder_parent + 1;
    
    SEXP r_postorder = PROTECT(Rf_allocVector(INTSXP, k));
    Memcpy(INTEGER(r_postorder), postorder, k);

    UNPROTECT(1);
    return r_postorder;
}


SEXP C_ephylo_preorder_tips(
    SEXP r_num_tips,
    SEXP r_num_nodes,
    SEXP r_root,
    SEXP r_right_child, 
    SEXP r_left_sib)
{
    int num_tips = *INTEGER(r_num_tips);
    int num_nodes = *INTEGER(r_num_nodes);
    int root = *INTEGER(r_root) - 1;
    int *right_child = INTEGER(r_right_child);
    int *left_sib = INTEGER(r_left_sib);

    int *stack = (int *) R_alloc(num_nodes, sizeof(*stack));
    
    int u;
    int v;
    int stack_top;

    int *preorder = (int *) R_alloc(num_nodes, sizeof(*preorder));

    stack_top = -1;
    for (u = (right_child[root]-1); u != -1; u = (left_sib[u]-1))
    {
        ++stack_top;
        stack[stack_top] = u;
    }

    int k = 0;
    while (stack_top >= 0)
    {
        u = stack[stack_top];
        if (u <= num_tips)
            preorder[k++] = u + 1;
        --stack_top;
        for (v = (right_child[u]-1); v != -1; v = (left_sib[v]-1))
        {
            ++stack_top;
            stack[stack_top] = v;
        }
    }
    
    SEXP r_preorder = PROTECT(Rf_allocVector(INTSXP, k));
    Memcpy(INTEGER(r_preorder), preorder, k);

    UNPROTECT(1);
    return r_preorder;
}
