#include <stdio.h>

/*= remove excess space in names */
void rmexspce(int numsp, char names[][11]);

/*= number of digits before decimal */
int ndig_intpart(double f);

/*= convert utreec to Newick format, possibly with internal edge labels */
/* NOTE: names changed: `excess' space removed from names */
char *ut2in_rf(int numsp, double *utreec, char *names, int *pn, char **ilabels);

/*= convert utreec (new) to Newick format */
void ut2n_rf(int numsp, double *utreec, char *names, int *pn, char *tstr);

/*= convert utreec (new) to Newick format, topology only */
void ut2nt_rf(int numsp, double *utreec, char *names, int *pn, char *tstr);

/*=pr_tstr_labels: newick tree with utreec internal edges labeled */
void pr_tstr_labels(FILE *otreefile, int ntaxa, double *utreec,
                    char* names, int *pn);