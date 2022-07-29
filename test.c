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
            int *pn, char (*inames)[11], int is_aa) {

  int i, j;

  if (type==1)
    fprintf(fp, "%d %d\n", ntaxa, nsites);

  for (j=0; j<ntaxa; j++) {
    if (type==0)
      fprintf(fp, ">%.*s (%s)\n", pn[j+1]-pn[j], names+pn[j], inames[j]);
    else
      fprintf(fp, "%.*s (%s) ", pn[j + 1] - pn[j], names + pn[j], inames[j]);
    for (i=0; i<nsites; i++) {
      if (is_aa) {
        fprintf(fp, "%c", i2lp(seq[i+j*nsites]));
      } else {
        fprintf(fp, "%c", i2l(seq[i+j*nsites]));
      }
    }
    fprintf(fp, "\n");
  }

  return;
}

void cleanup(FILE *fp, FILE *fpout, int *seq, char *names, int *pn, \
             char (*inames)[11]) {

  if (fp)     { fclose(fp);    fp=0;     }
  if (fpout)  { fclose(fpout); fpout=0;  }
  if (seq)    { free(seq);     seq=0;    }
  if (names)  { free(names);   names=0;  }
  if (pn)     { free(pn);      pn=0;     }
  if (inames) { free(inames);  inames=0; }
}

/*
 * test 1: read sequence from seq_fn (type 0 for fasta, 1 for phylip)
 * test 2: read Newick tree and replace taxon ids with digits
 */
void test_rseq(const char *seq_fn, const char *tree_fn, const int is_aa, const int type) {
  FILE *fp, *fpout;
  int *seq, ntaxa, nsites, *pn;
  char *names, (*inames)[11], *tree, *tmp_fn;

  fp = 0; fpout = 0; seq = 0; names = 0; tree = 0;

  if( !(type == 0 || type == 1) ) {
    printf("The type argument must either be 0 (fasta) or 1 (phylip)\n");
    return;
  }

  tmp_fn = malloc(strlen("./sample_data/") + strlen(seq_fn) + 1);
  strcpy(tmp_fn, "./sample_data/"); strcat(tmp_fn, seq_fn);
  if ( !(fp = fopen(tmp_fn, "r")) )
    printf("Failed to open the sequence file, %s for reading.\n", tmp_fn);
  free(tmp_fn); tmp_fn=0;
  tmp_fn = malloc(strlen("./out/") + strlen(seq_fn) + 1);
  strcpy(tmp_fn, "./out/"); strcat(tmp_fn, seq_fn);
  if ( !(fpout = fopen(tmp_fn, "w")) )
    printf("Failed to open file, %s for writing.\n", tmp_fn);
  free(tmp_fn); tmp_fn=0;
  if (type == 0) rseq_fasta(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames, is_aa);
  if (type == 1) rseq_rphy(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames, is_aa);
  palign(fpout, type, ntaxa, nsites, seq, names, pn, inames, is_aa);
  fclose(fp); fp=0; fclose(fpout); fpout=0;

  /* test converting tree to one with digits for taxon ids */
  if (tree_fn) {
    tmp_fn = malloc(strlen("./sample_data/") + strlen(tree_fn) + 1);
    strcpy(tmp_fn, "./sample_data/"); strcat(tmp_fn, tree_fn);
    if ( !(fp = fopen(tmp_fn, "r")) )
      printf("Failed to open the tree file, %s, for reading.\n", tmp_fn);
    free(tmp_fn); tmp_fn=0;
    itree(fp, &tree, names, pn, ntaxa);

    tmp_fn = malloc(strlen("./out/") + strlen(tree_fn) + 1);
    strcpy(tmp_fn, "./out/"); strcat(tmp_fn, tree_fn);
    if ( !(fpout = fopen(tmp_fn, "w")) )
      printf("Failed to open the itree output file, %s.\n", tmp_fn);
    free(tmp_fn); tmp_fn=0;
    if(tree) fprintf(fpout, "%s\n", tree);
  }
  cleanup(fp, fpout, seq, names, pn, inames);
  if (tree) { free(tree); tree=0; }

}

/* Read in a relaxed phylip sequence stored in seq_fn and a tree in Ed's utreec
 * format and print out the newick tree with and without branch lengths.
 */
void test_utreec(const char *seq_fn, const int is_aa, const char *utreec_fn) {
  char *names, (*inames)[11], *tree, tstr[TMAX], *tmp, *tmp_fn;
  double *utreec;
  FILE *fp, *fpout;
  int i, *seq, ntaxa, nsites, *pn;

  fp = 0; fpout = 0; seq = 0; names = 0; tree = 0;

  tmp_fn = malloc(strlen("./sample_data/") + strlen(seq_fn) + 1);
  strcpy(tmp_fn, "./sample_data/"); strcat(tmp_fn, seq_fn);
  if ( !(fp = fopen(tmp_fn, "r")) )
    printf("Failed to open the phylip sequence file, %s.\n", tmp_fn);
  free(tmp_fn); tmp_fn=0;

  rseq_rphy(fp, &ntaxa, &nsites, &seq, &names, &pn, &inames, is_aa);
  fclose(fp); fp=0;

  utreec = malloc((ntaxa*4)*sizeof(*utreec));
  tmp_fn = malloc(strlen("./sample_data/") + strlen(utreec_fn) + 1);
  strcpy(tmp_fn, "./sample_data/"); strcat(tmp_fn, utreec_fn);
  if ( !(fp = fopen(tmp_fn, "r")) )
    printf("Failed to open the utreec file, %s.\n", tmp_fn);
  free(tmp_fn); tmp_fn=0;
  for (i=0; i<(ntaxa-1)*4; i++)
    if ( fscanf(fp, "%lf", &utreec[i]) != 1 )
      printf("Failed to read a double from the utreec file.\n");
  fclose(fp); fp=0;

  ut2n_rf(ntaxa, utreec, names, pn, tstr);

  tmp_fn = malloc(strlen("./out/") + strlen(utreec_fn) + 7 + 1);
  strcpy(tmp_fn, "./out/"); strcat(tmp_fn, utreec_fn);
  tmp = strrchr(tmp_fn, '.'); if (!tmp) *tmp = '\0';
  strcat(tmp_fn, ".newick");
  if ( !(fpout = fopen(tmp_fn, "w")) )
    printf("Failed to open %s for writting.\n", tmp_fn);
  free(tmp_fn); tmp_fn=0;
  fprintf(fpout, "%s\n", tstr);
  fclose(fpout); fpout=0;

  ut2nt_rf(ntaxa, utreec, names, pn, tstr);

  tmp_fn = malloc(strlen("./out/") + strlen(utreec_fn) + 9 + 1);
  strcpy(tmp_fn, "./out/"); strcat(tmp_fn, utreec_fn);
  tmp = strrchr(tmp_fn, '.'); if (!tmp) *tmp = '\0';
  strcat(tmp_fn, "_t.newick");
  if ( !(fpout = fopen(tmp_fn, "w")) )
    printf("Failed to open %s for writting.\n", tmp_fn);
  free(tmp_fn); tmp_fn=0;
  fprintf(fpout, "%s\n", tstr);

  cleanup(fp, fpout, seq, names, pn, inames);
  if (tree) { free(tree); tree=0; }
  if (utreec) { free(utreec); utreec=0; }
}

int main() {

  /* With no tree file specified, just print the sequence file, but also with
     numeric taxon ids in parentheses.  The file is written under ./out. */
  test_rseq("1052.rd5_h0.7.bmge.fasta",0,1,0);

  /* This time with a relaxed phylip sequence and the associated tree.  The
     itree will also be written under ./out/. */
  test_rseq("1052.rd5_h0.7.bmge.phy","1052.rd5_lgc20.treefile",1,1);

  /* With a relaxed phylip sequence and the associated tree in Ed's utreec
     format, convert to the newick tree with and without branch lengths. Output
     is written to ./out/. */
  test_utreec("pmsf.phy", 1, "amborella.utreec");
  test_utreec("1052.rd5_h0.7.bmge.phy", 1, "1052.rd5_h0.7.bmge.utreec");

  /* test nucleotide sequences */
  test_rseq("lysin.nuc",0,0,1);
  test_rseq("lysin.nuc","lysin.tre",0,1);

  return(0);
}
