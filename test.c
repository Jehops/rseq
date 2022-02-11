#include <stdio.h>
#include "rseq.h"

int main() {

  //  int phe = l2ip('F');

  //printf("The value returned from l2ip() is: %d.\n", phe);

  FILE *fp;
  int *seq, ntaxa, nsites, *pn;
  char **names;

  fp = fopen("./1052.rd5_h0.7.bmge.fasta", "r");

  rseq_fasta(fp, &ntaxa, &nsites, seq=0, names=0, pn=0);

  printf("Number of taxa: %d\n", ntaxa);
  printf("Number of sites: %d\n", nsites);
  printf("Names:\n");
  for (int i=0; i<ntaxa; i++) {
    for (int j=0; j<pn[i]; j++) {
      printf("%c",names[i][j]);
    }
    printf("\n");
  }

  return(0);

}