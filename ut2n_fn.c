#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*= remove excess space in names */
void rmexspce(int numsp, char names[][11]) {
  int i,j;
  for(i = 0; i < numsp; i++){
    for(j = 9; j >= 0; j--)
      if (names[i][j] != ' ' && names[i][j] != '\0' &&  names[i][j] != '\t' && names[i][j] != '\r'){
	names[i][j+1] = '\0';
	break;
      }}
}

/*= number of digits before decimal */
int ndig_intpart(double f) {
  int d=1,ndigits=0;
  while(((int) f/d) > 0){
    d *= 10;
    ndigits++;
  }
  if (ndigits == 0) ndigits=1;
  return(ndigits);
}

/*= convert utreec to Newick format, possibly with internal edge labels */
/* NOTE: names changed: `excess' space removed from names */
char *ut2in_rf(int numsp, double *utreec, char *names, int *pn, char **ilabels) {

  char **tstring,*tstr;
  double x;
  int i, len, r, l;
  /* int alen; Only when alen in front of all sprintf() calls */

  tstring = (char **) malloc((2*numsp-2)*sizeof(*tstring));

  for(i = 0; i < numsp; i++) {
    len = pn[i+1]-pn[i];
    tstring[i] = malloc((len+1)*sizeof(**tstring));
    strncpy(tstring[i], names+pn[i], len);
    tstring[i][len] = '\0';
  }

  for(i = 0; i < (numsp-3); i++) {
    /* memory allocation */
    r = utreec[i*4]; l = utreec[1+i*4];
    len = strlen(tstring[r]) + strlen(tstring[l]) +
      ndig_intpart(utreec[2+i*4]) + ndig_intpart(utreec[3+i*4]) + 18;
    if (r >= numsp) len += strlen(ilabels[r-numsp]);
    if (l >= numsp) len += strlen(ilabels[l-numsp]);
    tstring[i+numsp] = malloc((len+1)*sizeof(**tstring));
    /* create subtree */
    if (r >= numsp && l >= numsp)
      sprintf(tstring[i+numsp],"(%s%s:%.5f,%s%s:%.5f)",tstring[r],
              ilabels[r-numsp],utreec[2+i*4],tstring[l],ilabels[l-numsp],
              utreec[3+i*4]);
    if (r >= numsp && l < numsp)
      sprintf(tstring[i+numsp],"(%s%s:%.5f,%s:%.5f)",tstring[r],
              ilabels[r-numsp],utreec[2+i*4],tstring[l],utreec[3+i*4]);
    if (r < numsp && l >= numsp)
      sprintf(tstring[i+numsp],"(%s:%.5f,%s%s:%.5f)",tstring[r],
              utreec[2+i*4],tstring[l],ilabels[l-numsp],utreec[3+i*4]);
    if (r < numsp && l < numsp)
      sprintf(tstring[i+numsp],"(%s:%.5f,%s:%.5f)",tstring[r],
              utreec[2+i*4],tstring[l],utreec[3+i*4]);
    /* printf("%i %i: %s\n",alen+1,i+numsp,tstring[i+numsp]); */
  }
  
  /* PHYLIP requires trifurcation */
  /* memory allocation */
  r = utreec[(numsp-3)*4]; l = utreec[1+(numsp-3)*4];
  len = strlen(tstring[r]) + strlen(tstring[l]) + 16 +
    ndig_intpart(utreec[2+(numsp-3)*4])+ndig_intpart(utreec[3+(numsp-3)*4]);
  if (r >= numsp) len += strlen(ilabels[r-numsp]);
  if (l >= numsp) len += strlen(ilabels[l-numsp]);
  tstring[2*numsp-3] = malloc((len+1)*sizeof(**tstring));
  /* create subtree */
  i=numsp-3;
  if (r >= numsp && l >= numsp)
    sprintf(tstring[i+numsp],"%s%s:%.5f,%s%s:%.5f",tstring[r],
            ilabels[r-numsp],utreec[2+i*4],tstring[l],ilabels[l-numsp],
            utreec[3+i*4]);
  if (r >= numsp && l < numsp)
    sprintf(tstring[i+numsp],"%s%s:%.5f,%s:%.5f",tstring[r],
            ilabels[r-numsp],utreec[2+i*4],tstring[l],utreec[3+i*4]);
  if (r < numsp && l >= numsp)
    sprintf(tstring[i+numsp],"%s:%.5f,%s%s:%.5f",tstring[r],
            utreec[2+i*4],tstring[l],ilabels[l-numsp],utreec[3+i*4]);
  if (r < numsp && l < numsp)
    sprintf(tstring[i+numsp],"%s:%.5f,%s:%.5f",tstring[r],
            utreec[2+i*4],tstring[l],utreec[3+i*4]);
  /* printf("%i %i %i %i %i: %s\n",r,l,len,alen+1,2*numsp-3,tstring[2*numsp-3]); */
  /* memory allocation for final string */
  r = utreec[(numsp-2)*4]; l = utreec[1+(numsp-2)*4];
  r = (r < l)?r:l;
  x=utreec[2+(numsp-2)*4]+utreec[3+(numsp-2)*4];
  len = strlen(tstring[2*numsp-3]) + strlen(tstring[r]) + 12 + 
    ndig_intpart(x);
  if (r >= numsp) len += strlen(ilabels[r-numsp]);
  tstr = malloc((len+1)*sizeof(*tstr));
  /* tree */
  if (r >= numsp)
    sprintf(tstr,"(%s,%s%s:%.5f);",tstring[2*numsp-3],tstring[r],
            ilabels[r-numsp],x);
  if (r < numsp)
    sprintf(tstr,"(%s,%s:%.5f);",tstring[2*numsp-3],tstring[r],x);
  tstr[len-1] = '\0';
  /* printf("%i %i %i %s\n",r,len,alen+1,tstr); */
  for(i = 0; i < 2*numsp-2; i++) free(tstring[i]);

  free(tstring);
  return(tstr);
}

/*= convert utreec (new) to Newick format */
void ut2n_rf(int numsp, double *utreec, char *names, int *pn, char *tstr) {

  char **tstring;
  int i,len,r,l;

  tstring = malloc((2*numsp-2)*sizeof(*tstring));

  for(i = 0; i < numsp; i++) {
    len = pn[i+1]-pn[i];
    tstring[i] = malloc((len+1)*sizeof(**tstring));
    strncpy(tstring[i], names+pn[i], len);
    tstring[i][len] = '\0';
  }

  for(i = 0; i < (numsp-3); i++) {
    r = utreec[i*4]; l = utreec[1+i*4];
    len = strlen(tstring[r]) + strlen(tstring[l])+23; /* space for branch lengths in the 100s */
    tstring[i+numsp] = malloc((len+1)*sizeof(**tstring));
    sprintf(tstring[i+numsp],"(%s:%.5f,%s:%.5f)",tstring[r],utreec[2+i*4],tstring[l],
	    utreec[3+i*4]);
  }
  
  /* PHYLIP requires trifurcation */
  r = utreec[(numsp-3)*4]; l = utreec[1+(numsp-3)*4];
  len = strlen(tstring[r]) + strlen(tstring[l]) + 21;
  tstring[2*numsp-3] = malloc((len+1)*sizeof(**tstring));
  sprintf(tstring[2*numsp-3],"%s:%.5f,%s:%.5f",tstring[r],utreec[2+(numsp-3)*4],tstring[l],
          utreec[3+(numsp-3)*4]);
  r = utreec[(numsp-2)*4]; l = utreec[1+(numsp-2)*4];
  r = (r < l)?r:l;
  sprintf(tstr,"(%s,%s:%.5f);",tstring[2*numsp-3],tstring[r],
          utreec[2+(numsp-2)*4]+utreec[3+(numsp-2)*4]);
  len = strlen(tstr);
  tstr[len] = '\0';
  for(i = 0; i < 2*numsp-2; i++) free(tstring[i]);
  free(tstring);
}

/*= convert utreec (new) to Newick format, topology only */
void ut2nt_rf(int numsp, double *utreec, char *names, int *pn, char *tstr) {

  char **tstring;
  int i,len,r,l;

  tstring = malloc((2*numsp-2)*sizeof(*tstring));

  for(i = 0; i < numsp; i++) {
    len = pn[i+1]-pn[i];
    tstring[i] = malloc((len+1)*sizeof(**tstring));
    strncpy(tstring[i], names+pn[i], len);
    tstring[i][len] = '\0';
  }

  for(i = 0; i < (numsp-3); i++) {
    r = utreec[i*4]; l = utreec[1+i*4];
    len = strlen(tstring[r]) + strlen(tstring[l])+4; 
    tstring[i+numsp] = malloc((len+1)*sizeof(**tstring));
    sprintf(tstring[i+numsp],"(%s,%s)",tstring[r],tstring[l]);
  }
  
  /* PHYLIP requires trifurcation */
  r = utreec[(numsp-3)*4]; l = utreec[1+(numsp-3)*4];
  len = strlen(tstring[r]) + strlen(tstring[l]) +2;
  tstring[2*numsp-3] = malloc((len+1)*sizeof(**tstring));
  sprintf(tstring[2*numsp-3],"%s,%s",tstring[r],tstring[l]);
  r = utreec[(numsp-2)*4]; l = utreec[1+(numsp-2)*4];
  r = (r < l)?r:l;
  sprintf(tstr,"(%s,%s);",tstring[2*numsp-3],tstring[r]);
  len = strlen(tstr);
  tstr[len] = '\0';
  for(i = 0; i < 2*numsp-2; i++) free(tstring[i]);
  free(tstring);
}

/*=pr_tstr_labels: newick tree with utreec internal edges labeled */
void pr_tstr_labels(FILE *otreefile, int ntaxa, double *utreec, 
		    char *names, int *pn) {
  
  char **ilabels,*tstr;
  int i,elab,len;

  ilabels = malloc((ntaxa-2)*sizeof(*ilabels));
  for(i = 0; i < ntaxa-2; i++) {
    elab=i+ntaxa;
    len=1;
    while(elab>=exp(len*log(10))) len++;
    ilabels[i] = malloc((len+1)*sizeof(**ilabels));
    sprintf(ilabels[i],"%i",elab);
  }
  tstr=ut2in_rf(ntaxa,utreec,names,pn,ilabels);
  fprintf(otreefile,"%s\n", tstr);

  for(i=0; i<ntaxa-2; i++) free(ilabels[i]);
  free(ilabels); free(tstr);
}
