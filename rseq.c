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

char ip2l(int a) {
  switch (a) {
  case alanine:        return('A');
  case arginine:       return('R');
  case asparagine:     return('N');
  case aspartic:       return('D');
  case cysteine:       return('C');
  case glutamine:      return('Q');
  case glutamic:       return('E');
  case glycine:        return('G');
  case histidine:      return('H');
  case isoleucine:     return('I');
  case leucine:        return('L');
  case lysine:         return('K');
  case methionine:     return('M');
  case phenylalanine:  return('F');
  case proline:        return('P');
  case serine:         return('S');
  case threonine:      return('T');
  case tryptophan:     return('W');
  case tyrosine:       return('Y');
  case valine:         return('V');
  case 21:             return('-');
  case 22:             return('X');
  case 23:             return('?');
  }

  return('?');
}

void _fasta_prescan(FILE *seqfp, int *ntaxa, int *nsites, int *totnl) {

  char cbuf;

  *totnl=0;
  while( (cbuf = fgetc(seqfp)) != EOF ) {
    if ( cbuf == '>' ) {
      (*ntaxa)++;
      while( (cbuf = fgetc(seqfp)) != '\n' ) {
        if (cbuf != '\r')
          (*totnl)++;
      }
    }
    else if ( *ntaxa == 1 ) {
      if ( cbuf != '\n' && cbuf != '\r' ) {
        (*nsites)++;
      }
    }
  }

  rewind(seqfp);
}

/*
 * Read a fasta alignment from the file pointed to by seqfp and store it in seq.
 * Amino acid values are determined from l2ip().  To retrieve the amino acid
 * at site i of species j, use seq[i + j*nsites].
 *
 * names:
 * pn: pn[j] gives length+1 of name[j], the name for species j

 * Allocated:
 *   seq: nsites*ntaxa chars
 *   names:
 *   pn:
 *

 */
void rseq_fasta(FILE *seqfp, int *ntaxa, int *nsites, int **seq,
                char **names, int **pn) {
  char cbuf, lbuf[LBUFLEN];
  int curseq=0, seqi=0, npi=0;

  *ntaxa = *nsites = 0;

  _fasta_ntaxa_nsites(seqfp, ntaxa, nsites);

  *names = malloc(*ntaxa * sizeof(char));
  *seq = malloc((*nsites)*(*ntaxa) * sizeof(int));
  *pn = malloc(*ntaxa * sizeof(int));

  while( (cbuf = fgetc(seqfp)) != EOF ) {
    if ( cbuf == '>' ) { // read sequence identifier line

      if (seqi) { // the previous sequence is in lbuf
        for(int i=0; i<seqi; i++) {
          (*seq)[i + curseq*(*nsites)] = l2ip(lbuf[i]);
        }
        seqi=0;
        curseq++;
      }

      if ( !fgets(lbuf, LBUFLEN, seqfp) ) { // read seq id into lbuf
        printf("Error reading a sequence name.\n");
      }
      lbuf[strcspn(lbuf, "\r\n")] = 0;
      strcpy((*names)+npi, lbuf);
      npi+=strlen(lbuf);
      (*pn)[curseq] = npi;
    } else { // read sequence
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
  for(int i=0; i<seqi; i++) {
    (*seq)[i + curseq*(*nsites)] = l2ip(lbuf[i]);
  }

  return;
}

void rseq_rphy(FILE *seqfile, int *ntaxa, int *nsites, int **seq, char **names,
               int **pn) {
  return;

}
