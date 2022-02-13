#include <stdio.h>
#include <stdlib.h>
#include "rseq.h"

int main() {

  FILE *fp=0, *fpout=0;
  int *seq=0, ntaxa, nsites, *pn=0;
  char *names=0;

  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.fasta", "r")) ) {
    printf("FATAL: Failed to open the sequence file.\n");
    exit(1);
  }

  if ( !(fpout = fopen("./rseq.out", "w")) ) {
    printf("FATAL: Failed to open the output file.\n");
    exit(1);
  }  

  rseq_fasta(fp, &ntaxa, &nsites, &seq, &names, &pn);

  printf("Number of taxa: %d\n", ntaxa);
  printf("Number of sites: %d\n", nsites);

  for (int i=0, p=0; i<ntaxa; p=pn[i], i++) {
    fprintf(fpout, ">%.*s\n", pn[i]-p, names+p);
    for (int j=0; j<nsites; j++) {
      fprintf(fpout, "%c", ip2l(seq[j+i*nsites]));
    }
    fprintf(fpout, "\n");
  }

  free(names);
  free(seq);
  free(pn);

  fclose(fp);
  fclose(fpout);

  return(0);

}