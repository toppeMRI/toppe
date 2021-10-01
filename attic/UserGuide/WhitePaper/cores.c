/*
  This file is part of the TOPPE development environment for platform-independent MR pulse programming.

  TOPPE is free software: you can redistribute it and/or modify
  it under the terms of the GNU Library General Public License as published by
  the Free Software Foundation version 2.0 of the License.

  TOPPE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Library General Public License for more details.

  You should have received a copy of the GNU Library General Public License
  along with TOPPE. If not, see <http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html>.
 
  (c) 2016 The Regents of the University of Michigan
  Jon-Fredrik Nielsen, jfnielse@umich.edu
 
  $Id: cores.c,v 1.29 2016/05/20 17:43:08 jfnielse Exp $
 */

#include <string.h>
#include <math.h>

#ifndef HW_IO
#include <stdio.h>
#else
#include <stdioLib.h>
#endif

#include "jfn_globaldefs.h"   
#include "cores.h"
#include "jfn_rf_files2.h"   /* need readshort() etc */


/*
 * Read .coredef file
 */
int cores_getinfo(char *fname, corestruct* coredefinfo) {

	FILE *fid;
	short nchar,tmp;
	int i;
	char line[200];

	fprintf(stderr,"cores_getinfo(): reading %s \n", fname);

	if ((fid = fopen (fname, "r")) == NULL)
		return(JFN_FAILURE);

	fgets(&line, 200, fid);   /* skip line */
	fscanf(fid, "%d\n", &(coredefinfo->ncores));
	fgets(&line, 200, fid);   /* skip line */

	coredefinfo->wavfiles = (char**)AllocNode(sizeof(char*)*(coredefinfo->ncores));
	coredefinfo->dur    = (int*)AllocNode(sizeof(int)*(coredefinfo->ncores));
	coredefinfo->hasRF  = (int*)AllocNode(sizeof(int)*(coredefinfo->ncores));
	coredefinfo->hasDAQ = (int*)AllocNode(sizeof(int)*(coredefinfo->ncores));

	for (i = 0; i < coredefinfo->ncores; i++)  {
		coredefinfo->wavfiles[i] = (char*)AllocNode(sizeof(char)*100);
		fscanf(fid, "%s\t%d\t%d\t%d\n", 
			coredefinfo->wavfiles[i], 
			&(coredefinfo->dur[i]),
			&(coredefinfo->hasRF[i]),
			&(coredefinfo->hasDAQ[i]));
		fprintf(stderr, "cores_getinfo(): %s\t%d\t%d\t%d\n", 
			coredefinfo->wavfiles[i], 
			coredefinfo->dur[i],
			coredefinfo->hasRF[i],
			coredefinfo->hasDAQ[i]);
	}

/*
	readshort(&(coredefinfo->ncores), 1, fid);
	fprintf(stderr, "\tcoredefinfo.ncores= %d \n", coredefinfo->ncores);

	for (i = 0; i < coredefinfo->ncores; i++)  {
		coredefinfo->wavfiles[i] = (char*)AllocNode(sizeof(char)*80);
		readshort(&nchar, 1, fid);
		fscanf(fid, "%s\n", coredefinfo->wavfiles[i]);
		readshort(&(coredefinfo->dur[i]), 1, fid);
		readshort(&(coredefinfo->hasRF[i]), 1, fid);
		readshort(&(coredefinfo->hasDAQ[i]), 1, fid);
	}
*/

	fclose (fid);
	return JFN_SUCCESS;
}


/* free memory */
int cores_freemem(corestruct *coredefinfo) {

	int i,c;
	int ncores = coredefinfo->ncores;

	FreeNode(coredefinfo->dur);
	FreeNode(coredefinfo->hasRF);
	FreeNode(coredefinfo->hasDAQ);

	for (i = 0; i < ncores; i++)  {
		FreeNode(coredefinfo->wavfiles[i]);
	}

	FreeNode(coredefinfo->wavfiles);

	return JFN_SUCCESS;
}

int cores_getloophdr(char *fname, int* loophdr) {

	FILE *fid;
	int nt, i;
	char line[200];

	fprintf(stderr,"cores_getloophdr(): reading %s \n", fname);

	if ((fid = fopen (fname, "r")) == NULL)
		return(JFN_FAILURE);

	fgets(&line, 200, fid);   /* go past first line */
	fscanf(fid, "%d\t%d\t%d\t%d\n", &(loophdr[0]), &(loophdr[1]), &(loophdr[2]), &(loophdr[3])); 
	fprintf(stderr, "\tcores_getloophdr: Found %d startseq() calls in loop file\n", loophdr[0]);

	fclose (fid);

	return JFN_SUCCESS;
}


/* read .loop file */
int cores_readloop(char *fname, int* looparr) {

	FILE *fid;
	int nt, i, j;
	char line[200];
	int loophdr[4];

	fprintf(stderr,"cores_readloop(): reading %s \n", fname);

	if ((fid = fopen (fname, "r")) == NULL)
		return(JFN_FAILURE);

	fgets(&line, 200, fid);   /* skip line */
	fscanf(fid, "%d\t%d\t%d\t%d\n", &(loophdr[0]), &(loophdr[1]), &(loophdr[2]), &(loophdr[3])); 
	nt = loophdr[0];

	fgets(&line, 200, fid);   /* skip line */

	/* load loop array */
	for (i = 0; i < nt; i++)  {
		for (j = 0; j < NL-1; j++) {	
			fscanf(fid, "%d\t", &(looparr[i*NL+j]));
		}
		fscanf(fid, "%d\n", &(looparr[i*NL+NL-1]));
/*
		fprintf(stderr,"cores_readloop(): phi = %d \n", looparr[i*NL+NL-1]);
		fscanf(fid, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			&(looparr[i*NL+0]), &(looparr[i*NL+1]), &(looparr[i*NL+2]), &(looparr[i*NL+3]), &(looparr[i*NL+4]), 
			&(looparr[i*NL+5]), &(looparr[i*NL+6]), &(looparr[i*NL+7]), &(looparr[i*NL+8]), &(looparr[i*NL+9]),
			&(looparr[i*NL+10]));
*/
	}

	for (i = nt-20; i < nt; i++)  {
		for (j = 0; j < NL-1; j++) {	
			fprintf(stderr, "%d\t", looparr[i*NL+j]);
		}
		fprintf(stderr, "%d\n", looparr[i*NL+NL-1]);

	}

	fclose (fid);

	return JFN_SUCCESS;
}

/* EOF */
