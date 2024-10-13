#include <R.h>
#include <Rinternals.h>
#include <assert.h>
#include <ctype.h>

/*
The Newick Tree Format grammar archived from:
https://phylipweb.github.io/phylip/newick_doc.html

Thursday, August 30, 1990

Gary Olsen's Interpretation of the "Newick's 8:45" Tree Format Standard

Conventions:
   Items in { } may appear zero or more times.
   Items in [ ] are optional, they may appear once or not at all.
   All other punctuation marks (colon, semicolon, parentheses, comma and
         single quote) are required parts of the format.


              tree ==> descendant_list [ root_label ] [ : branch_length ] ;

   descendant_list ==> ( subtree { , subtree } )

           subtree ==> descendant_list [internal_node_label] [: branch_length]
                   ==> leaf_label [: branch_length]

            root_label ==> label
   internal_node_label ==> label
            leaf_label ==> label

                 label ==> unquoted_label
                       ==> quoted_label

        unquoted_label ==> string_of_printing_characters
          quoted_label ==> ' string_of_printing_characters '

         branch_length ==> signed_number
                       ==> unsigned_number


Notes:
   Unquoted labels may not contain blanks, parentheses, square brackets,
        single_quotes, colons, semicolons, or commas.
   Underscore characters in unquoted labels are converted to blanks.
   Single quote characters in a quoted label are represented by two single
        quotes.
   Blanks or tabs may appear anywhere except within unquoted labels or
        branch_lengths.
   Newlines may appear anywhere except within labels or branch_lengths.
   Comments are enclosed in square brackets and may appear anywhere
        newlines are permitted.


Other notes:
   PAUP (David Swofford) allows nesting of comments.
   TreeAlign (Jotun Hein) writes a root node branch length (with a value of
        0.0).
   PHYLIP (Joseph Felsenstein) requires that an unrooted tree begin with a
        trifurcation; it will not "uproot" a rooted tree.


Example:

   (((One:0.2,Two:0.3):0.3,(Three:0.5,Four:0.3):0.2):0.3,Five:0.7):0.0;

           +-+ One
        +--+
        |  +--+ Two
     +--+
     |  | +----+ Three
     |  +-+
     |    +--+ Four
     +
     +------+ Five
*/

typedef struct newick_reader {
    int *parent;
    int *left_child;
    int *right_child;
    int *left_sib;
    int *right_sib;
    double *edge_len;
    SEXP node_label;
    PROTECT_INDEX node_label_idx;
    int bufsize;
    int curnode;
    int num_nodes;
    double brlen;
    char label[1024];
    FILE *newickfile;
} newick_reader_t;


static void
newick_free(newick_reader_t *reader)
{
    free(reader->parent);
    free(reader->left_child);
    free(reader->right_child);
    free(reader->left_sib);
    free(reader->right_sib);
    free(reader->edge_len);
    fclose(reader->newickfile);
    UNPROTECT(1); // reader->node_label
}


static void
newick_init(newick_reader_t *reader, const char *filename, int bufsize)
{
    reader->newickfile = fopen(filename, "r");
    if (!reader->newickfile)
        Rf_error("cannot open newick file");
    reader->bufsize = bufsize;
    reader->num_nodes = 0;
    reader->parent = malloc(bufsize * sizeof(int));
    reader->left_child = malloc(bufsize * sizeof(int));
    reader->right_child = malloc(bufsize * sizeof(int));
    reader->left_sib = malloc(bufsize * sizeof(int));
    reader->right_sib = malloc(bufsize * sizeof(int));
    reader->edge_len = malloc(bufsize * sizeof(double));
    PROTECT_WITH_INDEX(reader->node_label = Rf_allocVector(STRSXP, bufsize), 
        &reader->node_label_idx);
    if (!reader->parent || !reader->left_child || ! reader->right_child ||
        !reader->left_sib || !reader->right_sib || !reader->edge_len)
    {
        newick_free(reader);
        Rf_error("memory allocation failed");
    }
    memset(reader->label, 0, 1024);
}


static void
newick_maybe_realloc(newick_reader_t *reader)
{
    if (reader->num_nodes >= reader->bufsize)
    {
        reader->bufsize *= 2;
        int *parent = realloc(reader->parent, 
            sizeof(int) * reader->bufsize);
        int *left_child = realloc(reader->left_child, 
            sizeof(int) * reader->bufsize);
        int *right_child = realloc(reader->right_child, 
            sizeof(int) * reader->bufsize);
        int *left_sib = realloc(reader->left_sib, 
            sizeof(int) * reader->bufsize);
        int *right_sib = realloc(reader->right_sib, 
            sizeof(int) * reader->bufsize);
        double *edge_len = realloc(reader->edge_len, 
            sizeof(double) * reader->bufsize);
        if (!parent || !left_child || !right_child || !left_sib || 
            !right_sib || !edge_len)
        {
            newick_free(reader);
            Rf_error("memory allocation failed");
        }
        reader->parent = parent;
        reader->left_child = left_child;
        reader->right_child = right_child;
        reader->left_sib = left_sib;
        reader->right_sib = right_sib;
        reader->edge_len = edge_len;
        SEXP node_label = PROTECT(Rf_allocVector(STRSXP, reader->bufsize));
        for (int i = 0; i < reader->num_nodes; ++i)
        {
            SET_STRING_ELT(node_label, i, STRING_ELT(reader->node_label, i));
        }
        REPROTECT(node_label, reader->node_label_idx);
        reader->node_label = node_label;
        UNPROTECT(1);
    }
}


static char
newick_getc(newick_reader_t *reader)
{
    int c = fgetc(reader->newickfile);
    while (isspace(c))
        c = fgetc(reader->newickfile);
    return c;
}


static void
newick_ungetc(newick_reader_t *reader, char c)
{
    ungetc(c, reader->newickfile);
}


static void
newick_read_brlen(newick_reader_t *reader)
{
    char c = newick_getc(reader);
    if (c == ':')
    {
        int has_brlen = fscanf(reader->newickfile, "%lf", &reader->brlen);
        if (has_brlen)
            reader->edge_len[reader->curnode] = reader->brlen;
    }
    else
    {
        ungetc(c, reader->newickfile);
    }
}


static void
newick_read_label(newick_reader_t *reader)
{
    char c = newick_getc(reader);
    if (c == '[')
    {
        ungetc(c, reader->newickfile);
    }
    else
    {
        ungetc(c, reader->newickfile);
        int has_label = fscanf(
            reader->newickfile, " %1023[^,;:()]", reader->label);
        if (has_label)
        {
            int n = strlen(reader->label);
            while (isspace(reader->label[n-1]))
                --n;
            reader->label[n] = '\0';
            SET_STRING_ELT(reader->node_label, reader->curnode, 
                Rf_mkChar(reader->label));
        }
    }
}


static void
newick_read_comment(newick_reader_t *reader)
{
    char c = newick_getc(reader);
    if (c == '[')
    {
        fscanf(reader->newickfile, " %*[^]]");
        c = fgetc(reader->newickfile);
        if (c != ']')
            Rf_error("missing closing bracket in note");
    }
    else
    {
        ungetc(c, reader->newickfile);
    }
}


static void
newick_read(newick_reader_t *reader)
{
    char c;
    int u;
    int v;
    int paren = 0;
    fseek(reader->newickfile, -1, SEEK_END);
    c = fgetc(reader->newickfile);
    while (ftell(reader->newickfile) != SEEK_SET && isspace(c))
    {
        fseek(reader->newickfile, -2, SEEK_CUR);   
        c = fgetc(reader->newickfile);
    }
    if (c != ';')
    {
        newick_free(reader);
        Rf_error("missing terminating semicolon");
    }
    rewind(reader->newickfile);
    reader->curnode = 0;
    reader->num_nodes = 1;
    reader->left_child[0] = -1;
    reader->right_child[0] = -1;
    reader->left_sib[0] = -1;
    reader->right_sib[0] = -1;
    reader->parent[0] = -1;
    reader->edge_len[0] = NA_REAL;
    SET_STRING_ELT(reader->node_label, 0, NA_STRING);
    while ((c = newick_getc(reader)) != ';')
    {
        newick_maybe_realloc(reader);
        switch (c)
        {
            case '(': // starting a new clade
                u = reader->curnode;
                v = reader->num_nodes;
                reader->left_child[v] = -1;
                reader->right_child[v] = -1;
                reader->left_sib[v] = -1;
                reader->right_sib[v] = -1;
                reader->parent[v] = u;
                reader->edge_len[v] = NA_REAL;
                SET_STRING_ELT(reader->node_label, v, NA_STRING);
                if (reader->right_child[u] == -1)
                {
                    reader->left_child[u] = v;
                }
                else
                {
                    reader->left_sib[v] = reader->right_child[u];
                    reader->right_sib[reader->right_child[u]] = v;     
                }
                reader->right_child[u] = v;
                reader->curnode = v;
                reader->num_nodes += 1;
                ++paren;
                break;
            case ',': // starting a new clade sibling to the current clade
                u = reader->curnode;
                v = reader->num_nodes;
                reader->left_child[v] = -1;
                reader->right_child[v] = -1;
                reader->left_sib[v] = -1;
                reader->right_sib[v] = -1;
                reader->parent[v] = reader->parent[u];
                reader->edge_len[v] = NA_REAL;
                SET_STRING_ELT(reader->node_label, v, NA_STRING);
                if (reader->parent[u] == -1 ||
                    reader->left_child[reader->parent[u]] == -1 ||
                    reader->right_child[reader->parent[u]] != u)
                {
                    newick_free(reader);
                    Rf_error("bad newick string");
                }
                reader->left_sib[v] = u;
                reader->right_sib[u] = v;
                reader->right_child[reader->parent[u]] = v;
                reader->curnode = v;
                reader->num_nodes += 1;
                break;
            case ')': // finishing a clade
                u = reader->curnode;
                reader->curnode = reader->parent[u];
                if (reader->curnode == -1 ||
                    reader->left_child[reader->curnode] == -1 ||
                    reader->right_child[reader->curnode] == -1)
                {
                    newick_free(reader);
                    Rf_error("bad newick string");
                }
                --paren;
                break;
            default:
                newick_ungetc(reader, c);
                newick_read_label(reader);
                newick_read_comment(reader);
                newick_read_brlen(reader);
                break;
        }
    }
    if (paren != 0)
    {
        newick_free(reader);
        Rf_error("bad newick string");
    }
}

#define ISTIP(u) ((left_child[(u)] == -1 && right_child[(u)] == -1) ? 1 : 0)

SEXP
C_ephylo_read_newick(SEXP filename, SEXP bufsize)
{
    newick_reader_t reader;
    newick_init(&reader, CHAR(STRING_ELT(filename, 0)), *INTEGER(bufsize));
    newick_read(&reader);
    int num_nodes = reader.num_nodes;
    int stack_top;
    int *stack = (int *)R_alloc(num_nodes, sizeof(*stack));
    int *index_map = (int *)R_alloc(num_nodes, sizeof(*index_map));
    int *parent = reader.parent;
    int *left_child = reader.left_child;
    int *right_child = reader.right_child;
    int *left_sib = reader.left_sib;
    int *right_sib = reader.right_sib;
    double *edge_len = reader.edge_len;
    SEXP node_label = PROTECT(reader.node_label);
    int k;
    int u;
    int v;
    int r_index;
    int ntip = 0;
    int nnode = 0;
    stack_top = 0;
    stack[stack_top] = 0;
    // re-number node indices to match ape convention
    while (stack_top != -1)
    {
        u = stack[stack_top];
        --stack_top;
        for (v = right_child[u]; v != -1; v = left_sib[v])
        {
            ++stack_top;
            stack[stack_top] = v;
        }
        if (ISTIP(u))
            index_map[u] = ntip++;
        else
            index_map[u] = nnode++;
    }
    SEXP r_parent = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_left_child = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_right_child = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_left_sib = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_right_sib = PROTECT(Rf_allocVector(INTSXP, num_nodes));
    SEXP r_edge_len = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SEXP r_node_label = PROTECT(Rf_allocVector(STRSXP, num_nodes));
    SEXP r_ape_edge = PROTECT(Rf_allocMatrix(INTSXP, num_nodes-1, 2));
    SEXP r_time = PROTECT(Rf_allocVector(REALSXP, num_nodes));
    SET_REAL_ELT(r_time, ntip, 0);
    k = 0;
    stack_top = 0;
    stack[stack_top] = 0;
    while (stack_top != -1)
    {
        u = stack[stack_top];
        --stack_top;
        for (v = right_child[u]; v != -1; v = left_sib[v])
        {
            ++stack_top;
            stack[stack_top] = v;
        }
        if (ISTIP(u))
            r_index = index_map[u];
        else
            r_index = index_map[u] + ntip;
        if (parent[u] != -1)
            SET_INTEGER_ELT(r_parent, r_index, index_map[parent[u]] + 1 + ntip);
        else
            SET_INTEGER_ELT(r_parent, r_index, 0);
        if (left_child[u] != -1)
            SET_INTEGER_ELT(r_left_child, r_index, 
                index_map[left_child[u]] + 1 + (1 - ISTIP(left_child[u])) * ntip);
        else
            SET_INTEGER_ELT(r_left_child, r_index, 0);
        if (right_child[u] != -1)
            SET_INTEGER_ELT(r_right_child, r_index,
                index_map[right_child[u]] + 1 + (1 - ISTIP(right_child[u])) * ntip);
        else
            SET_INTEGER_ELT(r_right_child, r_index, 0);
        if (left_sib[u] != -1)
            SET_INTEGER_ELT(r_left_sib, r_index,
                index_map[left_sib[u]] + 1 + (1 - ISTIP(left_sib[u])) * ntip);
        else
            SET_INTEGER_ELT(r_left_sib, r_index, 0);
        if (right_sib[u] != -1)
            SET_INTEGER_ELT(r_right_sib, r_index,
                index_map[right_sib[u]] + 1 + (1 - ISTIP(right_sib[u])) * ntip);
        else
            SET_INTEGER_ELT(r_right_sib, r_index, 0);
        SET_REAL_ELT(r_edge_len, r_index, edge_len[u]);
        SET_STRING_ELT(r_node_label, r_index, STRING_ELT(node_label, u));
        if (parent[u] == -1) continue;
        SET_INTEGER_ELT(r_ape_edge, k, index_map[parent[u]] + 1 + ntip);
        SET_INTEGER_ELT(r_ape_edge, k + num_nodes - 1, r_index + 1);
        SET_REAL_ELT(
            r_time,
            r_index,
            REAL_ELT(r_time, index_map[parent[u]] + ntip) + edge_len[u]
        );
        ++k;
    }
    SEXP ans = PROTECT(Rf_allocVector(VECSXP, 11));
    SET_VECTOR_ELT(ans, 0, r_parent);
    SET_VECTOR_ELT(ans, 1, r_left_child);
    SET_VECTOR_ELT(ans, 2, r_right_child);
    SET_VECTOR_ELT(ans, 3, r_left_sib);
    SET_VECTOR_ELT(ans, 4, r_right_sib);
    SET_VECTOR_ELT(ans, 5, r_edge_len);
    SET_VECTOR_ELT(ans, 6, r_node_label);
    SET_VECTOR_ELT(ans, 7, r_ape_edge);
    SET_VECTOR_ELT(ans, 8, Rf_ScalarInteger(ntip));
    SET_VECTOR_ELT(ans, 9, Rf_ScalarInteger(nnode));
    SET_VECTOR_ELT(ans, 10, r_time);
    UNPROTECT(11);
    newick_free(&reader);
    return ans;
}

