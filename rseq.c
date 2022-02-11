#include <stdlib.h>
#include <string.h>

#include "rseq.h"

int l2ip(char c) {
  switch (c) {
  case 'A': return (alanine);
  case 'R': return (arginine);
  case 'N': return (asparagine);
  case 'D': return (aspartic);
  case 'C': return (cysteine);
  case 'Q': return (glutamine);
  case 'E': return (glutamic);
  case 'G': return (glycine);
  case 'H': return (histidine);
  case 'I': return (isoleucine);
  case 'L': return (leucine);
  case 'K': return (lysine);
  case 'M': return (methionine);
  case 'F': return (phenylalanine);
  case 'P': return (proline);
  case 'S': return (serine);
  case 'T': return (threonine);
  case 'W': return (tryptophan);
  case 'Y': return (tyrosine);
  case 'V': return (valine);
  case '-': return (21);
  case 'X': return (22);
  case '?': return (23);
  }

  return (100);
}

void _fasta_ntaxa_nsites(FILE *seqfp, int *ntaxa, int *nsites) {

  int nsitesp = 0;
  char cbuf;

  while( (cbuf = fgetc(seqfp)) != EOF ) {
    if ( cbuf == '>' ) {
      (*ntaxa)++;
      if ( *ntaxa > 1 && !nsitesp ) {
        nsitesp=1;
      }
      while( (cbuf = fgetc(seqfp)) != '\n' );
    }
    else if ( !(nsitesp) ) {
      if ( cbuf != '\n' && cbuf != '\r' ) {
        (*nsites)++;
        printf("%d: %c\n", *nsites,cbuf);
      }
    }
  }
}

/*
 * seq: Allocate nsites*ntaxa spaces based on seqfile: seq[i + j*nsites].
 *      The amino acid value for site i species j comes from l2ip().
 *
 * pn: Length ntaxa.  pnames[j] gives position of the last character for species j.
 */
void rseq_fasta(FILE *seqfp, int *ntaxa, int *nsites, int *seq,
                char **names, int *pn) {
  char cbuf;
  char lbuf[LBUFLEN];
  int curseq=0;
  int seqi=0;

  *ntaxa = *nsites = 0;

  _fasta_ntaxa_nsites(seqfp, ntaxa, nsites);
  names = malloc(*ntaxa*sizeof(char));
  seq = malloc((*nsites)*(*ntaxa)*sizeof(char));
  pn = malloc(*ntaxa*sizeof(int));

  rewind(seqfp);

  while( (cbuf = fgetc(seqfp)) != EOF ) {
    if ( cbuf == '>' ) {

      if (seqi) { // there is a sequence in the buffer
        for(int i=0; i<seqi; i++) {
          seq[curseq*(*nsites)+i] = lbuf[i];
        }
        seqi=0;
      }

      fgets(lbuf, LBUFLEN, seqfp);
      names[curseq] = malloc(strlen(lbuf)*sizeof(char));
      strncpy(names[curseq],lbuf,LBUFLEN);
      pn[curseq++] = strnlen(lbuf,LBUFLEN);
    } // cbuf == '>'
    else {
      if ( cbuf != '\n' && cbuf != '\r' ) {
        if (seqi < LBUFLEN) {
          lbuf[seqi++] = cbuf;
        }
        else {
          printf("Buffer overflow\n");
        }
      }
    }
  }

  return;
}

void rseq_rphy(FILE *seqfile, int *ntaxa, int *nsites, int *seq, char *names[],
               int *pn) {
  return;

}
