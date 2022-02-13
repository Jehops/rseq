#include <stdio.h>
#include <stdlib.h>
#include "rseq.h"

/*
 * print an alignment
 * type 0: fasta
 * type 1: relaxed phylip
 */
void palign(FILE *fp, int type, int ntaxa, int nsites, char *names,
            int *pn, int* seq) {
  if (type==1)
    fprintf(fp, "%d %d\n", ntaxa, nsites);

  for (int j=0; j<ntaxa; j++) {
    if (type==0)
      fprintf(fp, ">%.*s\n", pn[j+1]-pn[j], names+pn[j]);
    else
      fprintf(fp, "%.*s ", pn[j + 1] - pn[j], names + pn[j]);
    for (int i=0; i<nsites; i++)
      fprintf(fp, "%c", ip2l(seq[i+j*nsites]));
    fprintf(fp, "\n");
  }

  return;
}

void cleanup(char *names, int *seq, int *pn, FILE *fp, FILE *fpout) {
  free(names);
  free(seq);
  free(pn);

  fclose(fp);
  fclose(fpout);
}


int main() {

  FILE *fp=0, *fpout=0;
  int *seq=0, ntaxa, nsites, *pn=0;
  char *names=0;

  // test fasta
  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.fasta", "r")) ) {
    printf("FATAL: Failed to open the fasta sequence file.\n");
    exit(1);
  }
  if ( !(fpout = fopen("./rseq_fasta.out", "w")) ) {
    printf("FATAL: Failed to open the fasta output file.\n");
    exit(1);
  }

  rseq_fasta(fp, &ntaxa, &nsites, &seq, &names, &pn);
  palign(fpout, 0, ntaxa, nsites, names, pn, seq);

  cleanup(names, seq, pn, fp, fpout);

  // test relaxed phylip
  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.phy", "r")) ) {
    printf("FATAL: Failed to open the phylip sequence file.\n");
    exit(1);
  }
  if ( !(fpout = fopen("./rseq_phylip.out", "w")) ) {
    printf("FATAL: Failed to open the phylip output file.\n");
    exit(1);
  }

  rseq_rphy(fp, &ntaxa, &nsites, &seq, &names, &pn);
  palign(fpout, 1, ntaxa, nsites, names, pn, seq);

  cleanup(names, seq, pn, fp, fpout);

  return(0);
}
