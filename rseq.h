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

int l2ip(char c);

char ip2l(int a);

void _fasta_prescan(FILE *seqfp, int *ntaxa, int *nsites, int *totnl);

void rseq_fasta(FILE *seqfile, int *ntaxa, int *nsites, int **seq, char **names,
                int **pn);

void rseq_rphy(FILE *seqfile, int *ntaxa, int *nsites, int **seq, char **names,
               int **pn);
