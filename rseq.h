#include <stdio.h>

#define LBUFLEN 4096

enum {
  alanine = 0,
  arginine,
  asparagine,
  aspartic,
  cysteine,
  glutamine,
  glutamic,
  glycine,
  histidine,
  isoleucine,
  leucine,
  lysine,
  methionine,
  phenylalanine,
  proline,
  serine,
  threonine,
  tryptophan,
  tyrosine,
  valine
};

/* letter to integer */
int l2i(char c);

/* integer to letter */
char i2l(int i);

/* letter to integer protein */
int l2ip(char c);

/* integer to letter protein */
char i2lp(int a);

void _fasta_prescan(FILE *seqfp, int *ntaxa, int *nsites, int *totnl);

void rseq_fasta(FILE *seqfile, int *ntaxa, int *nsites, int **seq, char **names,
                int **pn, char (**inames)[11], int is_aa);

void rseq_rphy(FILE *seqfp, int *ntaxa, int *nsites, int **seq, char **names,
               int **pn, char (**inames)[11], int is_aa);

void itree(FILE *treefp, char **itree, const char *names, const int *pn,
           const int ntaxa);