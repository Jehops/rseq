#include <stdio.h>
#include <stdlib.h>
#include "rseq.h"

/*
 * print an alignment
 * type 0: fasta
 * type 1: relaxed phylip
 */
void palign(FILE *fp, int type, int ntaxa, int nsites, int* seq, char *names, \
            int *pn, char (*inames)[11]) {
  if (type==1)
    fprintf(fp, "%d %d\n", ntaxa, nsites);

  for (int j=0; j<ntaxa; j++) {
    if (type==0)
      fprintf(fp, ">%.*s (%s)\n", pn[j+1]-pn[j], names+pn[j], inames[j]);
    else
      fprintf(fp, "%.*s (%s) ", pn[j + 1] - pn[j], names + pn[j], inames[j]);
    for (int i=0; i<nsites; i++)
      fprintf(fp, "%c", ip2l(seq[i+j*nsites]));
    fprintf(fp, "\n");
  }

  return;
}

void cleanup(FILE *fp, FILE *fpout, int *seq, char *names, int *pn, \
             char (*inames)[11]) {

  fclose(fp);
  fclose(fpout);

  free(seq);
  free(names);
  free(pn);
  free(inames);
}


int main() {

  FILE *fp=0, *fpout=0;
  int *seq=0, ntaxa, nsites, *pn=0;
  char *names=0, (*inames)[11], *tree=0;

  // test fasta
  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.fasta", "r")) ) {
    printf("FATAL: Failed to open the fasta sequence file.\n");
    exit(1);
  }
  if ( !(fpout = fopen("./rseq_fasta.out", "w")) ) {
    printf("FATAL: Failed to open the fasta output file.\n");
    exit(1);
  }

  rseq_fasta(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames);
  palign(fpout, 0, ntaxa, nsites, seq, names, pn, inames);

  cleanup(fp, fpout, seq, names, pn, inames);

  // test relaxed phylip
  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.phy", "r")) ) {
    printf("FATAL: Failed to open the phylip sequence file.\n");
    exit(1);
  }
  if ( !(fpout = fopen("./rseq_phylip.out", "w")) ) {
    printf("FATAL: Failed to open the phylip output file.\n");
    exit(1);
  }

  rseq_rphy(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames);
  palign(fpout, 1, ntaxa, nsites, seq, names, pn, inames);

  fclose(fp);

  /* test converting tree to one with digits for taxon ids */
  if ( !(fp = fopen("./sample_data/1052.rd5_lgc20.treefile", "r")) ) {
    printf("FATAL: Failed to open the tree file.\n");
    exit(1);
  }
  itree(fp, &tree, names, pn, ntaxa);


  printf("%s\n",tree);

  cleanup(fp, fpout, seq, names, pn, inames);
  free(tree);
  fclose(fp);

  return(0);
}
