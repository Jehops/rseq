#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rseq.h"
#include "ut2n_fn.h"

#define TMAX 1048576 /* Account for tree files up to 1 MiB */

/*
 * print an alignment
 * type 0: fasta
 * type 1: relaxed phylip
 */
void palign(FILE *fp, int type, int ntaxa, int nsites, int* seq, char *names,
            int *pn, char (*inames)[11]) {

  int i, j;

  if (type==1)
    fprintf(fp, "%d %d\n", ntaxa, nsites);

  for (j=0; j<ntaxa; j++) {
    if (type==0)
      fprintf(fp, ">%.*s (%s)\n", pn[j+1]-pn[j], names+pn[j], inames[j]);
    else
      fprintf(fp, "%.*s (%s) ", pn[j + 1] - pn[j], names + pn[j], inames[j]);
    for (i=0; i<nsites; i++)
      fprintf(fp, "%c", i2lp(seq[i+j*nsites]));
    fprintf(fp, "\n");
  }

  return;
}

void cleanup(FILE *fp, FILE *fpout, int *seq, char *names, int *pn, \
             char (*inames)[11]) {

  if (fp) fclose(fp); if (fpout) fclose(fpout);
  if (seq) free(seq); if (names) free(names);
  if (pn) free(pn); if (inames) free(inames);
}

void test_rseq() {
  FILE *fp, *fpout;
  int *seq, ntaxa, nsites, *pn;
  char *names, (*inames)[11], *tree;

  fp = 0; fpout = 0; seq = 0; names = 0; tree = 0;

  /* test fasta */
  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.fasta", "r")) ) {
    printf("FATAL: Failed to open the fasta sequence file.\n");
    exit(1);
  }
  if ( !(fpout = fopen("./out/rseq_fasta.out", "w")) ) {
    printf("FATAL: Failed to open the fasta output file.\n");
    exit(1);
  }

  rseq_fasta(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames);
  palign(fpout, 0, ntaxa, nsites, seq, names, pn, inames);

  cleanup(fp, fpout, seq, names, pn, inames);

  /* test relaxed phylip */
  if ( !(fp = fopen("./sample_data/1052.rd5_h0.7.bmge.phy", "r")) ) {
    printf("FATAL: Failed to open the phylip sequence file.\n");
    exit(1);
  }
  if ( !(fpout = fopen("./out/rseq_phylip.out", "w")) ) {
    printf("FATAL: Failed to open the phylip output file.\n");
    exit(1);
  }

  rseq_rphy(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames);
  palign(fpout, 1, ntaxa, nsites, seq, names, pn, inames);

  fclose(fp); fclose(fpout);

  /* test converting tree to one with digits for taxon ids */
  if ( !(fp = fopen("./sample_data/1052.rd5_lgc20.treefile", "r")) ) {
    printf("FATAL: Failed to open the tree file.\n");
    exit(1);
  }
  itree(fp, &tree, names, pn, ntaxa);

  if ( !(fpout = fopen("./out/rseq_itree.out", "w")) ) {
    printf("FATAL: Failed to open the itree output file.\n");
    exit(1);
  }
  if(tree) fprintf(fpout, "%s\n", tree);

  cleanup(fp, fpout, seq, names, pn, inames);
  if (tree) free(tree);
}

void test_utree(const char *seq_fn, const char *utree_fn) {
  char *names, (*inames)[11], *tree, tstr[TMAX], *tmp, *tmp_fn;
  double *utreec;
  FILE *fp, *fpout;
  int i, *seq, ntaxa, nsites, *pn;

  fp = 0; fpout = 0; seq = 0; names = 0; tree = 0;

  /* Test amborella data */
  tmp_fn = malloc(strlen("./sample_data/") + strlen(seq_fn) + 1);
  strcpy(tmp_fn, "./sample_data/"); strcat(tmp_fn, seq_fn);
  if ( !(fp = fopen(tmp_fn, "r")) ) {
    printf("FATAL: Failed to open the phylip sequence file, %s.\n", tmp_fn);
    exit(1);
  }
  free(tmp_fn);

  rseq_rphy(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames);
  fclose(fp); fp=0;

  utreec = malloc((ntaxa*4)*sizeof(*utreec));
  tmp_fn = malloc(strlen("./sample_data/") + strlen(utree_fn) + 1);
  strcpy(tmp_fn, "./sample_data/"); strcat(tmp_fn, utree_fn);
  if ( !(fp = fopen(tmp_fn, "r")) ) {
    printf("FATAL: Failed to open the utree file, %s.\n", tmp_fn);
    exit(1);
  }
  free(tmp_fn);
  for (i=0; i<(ntaxa-1)*4; i++) {
    if ( fscanf(fp, "%lf", &utreec[i]) != 1 ) {
      printf("FATAL: Failed to read a double from the utreec file.\n");
      exit(1);
    }
  }
  fclose(fp); fp=0;

  ut2n_rf(ntaxa, utreec, names, pn, tstr);

  tmp_fn = malloc(strlen("./out/") + strlen(utree_fn) + 7 + 1);
  strcpy(tmp_fn, "./out/"); strcat(tmp_fn, utree_fn);
  tmp = strrchr(tmp_fn, '.'); if (!tmp) *tmp = '\0';
  strcat(tmp_fn, ".newick");
  if ( !(fpout = fopen(tmp_fn, "w")) ) {
    printf("FATAL: Failed to open %s for writting.\n", tmp_fn);
    exit(1);
  }
  free(tmp_fn);
  fprintf(fpout, "%s\n", tstr);
  fclose(fpout); fpout=0;

  ut2nt_rf(ntaxa, utreec, names, pn, tstr);

  tmp_fn = malloc(strlen("./out/") + strlen(utree_fn) + 9 + 1);
  strcpy(tmp_fn, "./out/"); strcat(tmp_fn, utree_fn);
  tmp = strrchr(tmp_fn, '.'); if (!tmp) *tmp = '\0';
  strcat(tmp_fn, "_t.newick");
  if ( !(fpout = fopen(tmp_fn, "w")) ) {
    printf("FATAL: Failed to open %s for writting.\n", tmp_fn);
    exit(1);
  }
  free(tmp_fn);
  fprintf(fpout, "%s\n", tstr);

  cleanup(fp, fpout, seq, names, pn, inames);
  if (tree) free(tree); if (utreec) free(utreec);
}

int main() {

  test_rseq();
  test_utree("pmsf.phy", "amborella.utreec");
  test_utree("1052.rd5_h0.7.bmge.phy", "1052.rd5_h0.7.bmge.utreec");

  return(0);
}
