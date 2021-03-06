#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "rseq.h"

int l2i(char c) {
  if (c == 'a' || c == 'A') return(0);
  if (c == 'c' || c == 'C') return(1);
  if (c == 'g' || c == 'G') return(2);
  if (c == 't' || c == 'T') return(3);
  if (c == '-')             return(4);

  return(5);
}

char i2l(int i) {
  switch(i) {
  case 0: return('A');
  case 1: return('C');
  case 2: return('G');
  case 3: return('T');
  case 4: return('-');
  }

  return('E');
}


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

char i2lp(int a) {
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

/* Read a fasta alignment of either nucleotides or amino acids from the file
 * pointed to by seqfp and store it in seq.  Nucleotide values are determined
 * from l2i() and amino acid values are determined from l2ip().  To retrieve the
 * amino acid at site i for taxon j, use seq[i + j*nsites].
 *
 * ntaxa:  number of taxa
 * nsites: number of aa sites
 * seq:    aa sequence data
 * names:  all taxon ids in a flat character array
 * pn:     pn[i] gives the starting index for ith taxon id
 * inames: integer taxon ids (string 10 characters long padded with blanks)
 * is_aa:  0 if working with nucleotides, 1 if working with amino acids
 *
 * Allocated with caller responsibility to free:
 *         seq: nsites*ntaxa chars
 *         names: total length all taxon ids chars
 *         pn: ntaxa+1 integers
 *         inames: ntaxa * (char*)[11]
 */
void rseq_fasta(FILE *seqfp, int *ntaxa, int *nsites, int **seq, char **names,
                int **pn, char (**inames)[11], int is_aa) {

  char cbuf, lbuf[LBUFLEN];
  int curseq=0, npi=0, seqi=0, i, totnl=0;

  *ntaxa = *nsites = 0;

  _fasta_prescan(seqfp, ntaxa, nsites, &totnl);

  *names = malloc((totnl+1)*sizeof(**names));
  *seq = malloc((*nsites)*(*ntaxa)*sizeof(**seq));
  *pn = malloc((*ntaxa+1)*sizeof(**pn));
  *inames = malloc(*ntaxa*sizeof((*inames)[11]));

  while( (cbuf = fgetc(seqfp)) != EOF ) {
    if ( cbuf == '>' ) { /* read sequence identifier line */

      if (seqi) { /* the previous sequence is in lbuf */
        for(i=0; i<seqi; i++) {
          if (is_aa) {
            (*seq)[i + curseq*(*nsites)] = l2ip(lbuf[i]);
          } else {
            (*seq)[i + curseq*(*nsites)] = l2i(lbuf[i]);
          }
        }
        seqi=0;
        curseq++;
      }

      if ( !fgets(lbuf, LBUFLEN, seqfp) ) { /* read seq id into lbuf */
        printf("Error reading a sequence name.\n");
      }
      lbuf[strcspn(lbuf, " 	\r\n")] = '\0';
      strcpy((*names)+npi, lbuf);
      (*pn)[curseq] = npi;
      npi+=strlen(lbuf);
      if ( sprintf(lbuf, "%d", curseq) > 9 )
        printf("%s:%d: Buffer overflow.\n", __FILE__, __LINE__-1);
      strcpy((*inames)[curseq],lbuf);
      for (i=strlen((*inames)[curseq]); i<10; i++)
        (*inames)[curseq][i] = ' ';
      (*inames)[curseq][10] = '\0';
    } else { /* read sequence */
      if ( cbuf != '\n' && cbuf != '\r' ) {
        if (seqi < LBUFLEN) {
          lbuf[seqi++] = cbuf;
        }
        else {
          printf("Buffer overflow.\n");
        }
      }
    }
  }
  for(i=0; i<seqi; i++) {
    if (is_aa) {
      (*seq)[i + curseq * (*nsites)] = l2ip(lbuf[i]);
    } else {
      (*seq)[i + curseq * (*nsites)] = l2i(lbuf[i]);
    }
  }
  (*pn)[++curseq] = npi;

  return;
}

/* Read a relaxed phylip alignment of either nucleotides or amino acids from the
 * file pointed to by seqfp and store it in seq.  Nucleotide values are
 * determined from l2i() and amino acid values are determined from l2ip().  To
 * retrieve the amino acid at site i for taxon j, use seq[i + j*nsites].
 *
 * ntaxa:  number of taxa
 * nsites: number of aa sites
 * seq:    nucleotide or aa sequence data
 * names:  all taxon ids in a flat character array
 * pn:     pn[i] gives the starting index for ith taxon id
 * inames: integer taxon ids (string 10 characters long padded with blanks)
 * is_aa:  0 if working with nucleotides, 1 if working with amino acids
 *
 * Allocated with caller responsibility to free:
 *         seq: nsites*ntaxa chars
 *         names: total length all taxon ids chars
 *         pn: ntaxa+1 integers
 *         inames: ntaxa * (char*)[11]
 */
void rseq_rphy(FILE *seqfp, int *ntaxa, int *nsites, int **seq, char **names,
               int **pn, char (**inames)[11], int is_aa) {

  char cbuf, lbuf[LBUFLEN];
  char **tnames = 0; /* temporary names array */
  int i, j, namel, npi=0, tnamel = 0; /* tnamel: total name length */

  if ( fscanf(seqfp, "%d %d", ntaxa, nsites) != 2 ) {
    printf("FATAL: Failed to read the phylip header.\n");
    return;
  }

  tnames = malloc((*ntaxa)*sizeof(char*));
  *seq = malloc((*nsites)*(*ntaxa)*sizeof(int));
  *pn = malloc((*ntaxa+1)*sizeof(int));
  *inames = malloc(*ntaxa*sizeof((*inames)[11]));

  for (j=0, namel=0; j<*ntaxa; j++) {
    if ( fscanf(seqfp, "%s", lbuf) != 1 ) { /* taxon id */
      printf("FATAL: Failed to read a taxon id.\n");
      free(tnames);  tnames = 0;
      free(*seq);    *seq = 0;
      free(*pn);     *pn = 0;
      free(*inames); *inames = 0;
      return;
    }
    namel = strlen(lbuf);
    tnames[j] = malloc((namel+1)*sizeof(char));
    tnamel+=namel;
    strcpy(tnames[j], lbuf);

    /* sequence */
    i=0;
    while( (cbuf = fgetc(seqfp)) != '\n' ) {
      if( ! isspace(cbuf) ) {
        if (is_aa) {
          (*seq)[i + j*(*nsites)] = l2ip(cbuf);
        } else {
          (*seq)[i + j*(*nsites)] = l2i(cbuf);
        }
        i++;
      }
    }
  }

  *names = malloc((tnamel+1)*sizeof(**names));
  for (j=0; j<*ntaxa; j++) {
    strcpy((*names)+npi, tnames[j]);
    (*pn)[j] = npi;
    npi+=strlen(tnames[j]);
    free(tnames[j]); tnames[j]=0;

    if ( sprintf(lbuf, "%d", j) > 9 )
      printf("%s:%d: Buffer overflow.\n", __FILE__, __LINE__-1);
    strcpy((*inames)[j],lbuf);
    for (i=strlen((*inames)[j]); i<10; i++)
      (*inames)[j][i] = ' ';
    (*inames)[j][10] = '\0';
  }
  (*pn)[*ntaxa] = npi;

  free(tnames); tnames=0;

  return;
}

/* Is char at tree[pos] after a '(' or ','? */
int afterpc(const char *tree, const int pos) {
  int i;

  if (pos == 0) return 1;

  for(i=pos-1; i>=0; i--) {
    if ( tree[i] == '(' || tree[i] == ',' )
      return 1;
    else if ( isspace(tree[i]) )
      ;
    else
      return 0;
  }
  return 1;
}

/* Is c a character that could be part of a taxon id? */
int idchar(char c) {
  if ( !isspace(c) && c != ':' && c != ';' && c != '(' && c != ')' && c != '['
       && c != ']' && c != ',' )
    return 1;
  return 0;
}

int idatpos(const int *pn, const int ntaxa, const int pos) {
  int i;

  for (i=0; i<ntaxa; i++) {
    if (pn[i] == pos)
      return i;
  }

  return -1;
}

/* Read a Newick tree from the file pointed to by treefp and store and replace
 * the taxon ids with ids that are integers determined from names and pn.
 *
 * itree: integer id tree
 * names: All taxon ids in a flat character array. If null, set numeric id as
 *        order encountered in tree.
 * pn:    pn[i] gives the starting index for ith taxon id
 * ntaxa: number of taxa
 *
 * Allocated with caller responsibility to free:
 *        itree
 */
void itree(FILE *treefp, char **itree, const char *names, const int *pn,
           const int ntaxa) {

  int ttl = 0; /* ttl is total tree length */
  int  i, id, j, k, l, pos, ridf; /* ridf is reading id flag */
  char cbuf, *ntree = 0, nbuff[LBUFLEN];

  /* read in the tree */
  while( (fgetc(treefp)) != EOF ) ttl++;
  ntree = malloc(ttl*sizeof(*ntree)+1);
  *itree = malloc(ttl*sizeof(**itree)+1);
  rewind(treefp);
  i=0;
  while( (cbuf = fgetc(treefp)) != EOF ) {
    ntree[i++] = cbuf;
  }
  ntree[i] = '\0';

  id=0; k=0; ridf=0;
  for (i=0; i<ttl; i++) {
    if ( idchar(ntree[i]) ) { /* ntree[i] could be part of an id */
      if (!ridf) {
        if (afterpc(ntree, i)) {
          ridf=1;
          j=0;
          nbuff[j++] = ntree[i];
        } else { /* This is not part of a taxon id */
          (*itree)[k++] = ntree[i];
        }
      } else {
        nbuff[j++] = ntree[i];
      }
    } else {
      if (ridf) { /* done reading a taxon id */
        nbuff[j++] = '\0';
        if (names != 0) { /* use names to get the numeric id */
          pos = strstr(names, nbuff) - names;
          if (pos < 0) {
            printf("%s:%d: There was an error retrieving the taxon name %s.\n",
                   __FILE__, __LINE__, nbuff);
            free(ntree);  ntree = 0;
            free(*itree); *itree = 0;
            return;
          }
          id = idatpos(pn, ntaxa, pos);
          l = sprintf((*itree)+k,"%d", id);
          if (l>9) printf("%s:%d: Buffer overflow.\n", __FILE__, __LINE__-1);
          k += l;
        } else { /* set numeric ids as order encountered in tree */
          l = sprintf((*itree)+k, "%d", id++);
          if (l>9) printf("%s:%d: Buffer overflow.\n", __FILE__, __LINE__-1);
          k += l;
        }
        ridf = 0;
      }
      (*itree)[k++] = ntree[i];
    }
  }

  (*itree)[k] = '\0';
  free(ntree); ntree=0;

  return;
}
