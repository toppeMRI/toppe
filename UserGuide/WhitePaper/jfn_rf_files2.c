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

  $Id: jfn_rf_files2.c,v 1.8 2017/08/29 17:08:37 jfnielse Exp $
 */

#include <string.h>
#include <math.h>

#ifndef HW_IO
#include <stdio.h>
#else
#include <stdioLib.h>
#endif

#include "jfn_globaldefs.h"   
#include "jfn_rf_files2.h"


/* get name of rf file from the names listed in 'rffiles.txt' */
int jfn_rf_getfilename(char *fname, int index, const char* directory) {
	char tmpfname[80];
	int i;
	FILE *fid;

  	sprintf(tmpfname, "%srffiles.txt", directory);
	if ((fid = fopen (tmpfname, "r")) == NULL)
		return(JFN_FAILURE);

	for (i = 0; i<index; i++)
		fscanf(fid, "%s\n", tmpfname);

	fscanf(fid, "%s\n", fname);

	fclose (fid);
	return JFN_SUCCESS;
}


short orderbyte(short in)
{
#ifdef LITTLE_ENDIAN

  char* sw;
  char tmp;
  sw = (char*)&in;
  tmp=sw[0];
  sw[0] = sw[1];
  sw[1] = tmp;

  /*
  in = htons(in);
  */
#endif
  return in;
}

unsigned short orderbyteu(unsigned short in)
{
#ifdef LITTLE_ENDIAN
  char* sw;
  char tmp;
  sw = (char*)&in;
  tmp=sw[0];
  sw[0] = sw[1];
  sw[1] = tmp;

  /*
  in = htons(in);
  */
#endif
  return in;
}


int jfn_rf_readheader(char *fname, rfstruct *rfinfo)
{
	FILE *fid;
	short nchar;
	int i;

	fprintf(stderr,"jfn_rf_readheader(): reading %s \n", fname);

	strcpy(rfinfo->fname, fname);

	if ((fid = fopen (fname, "r")) == NULL)
		return(JFN_FAILURE);

	/* go past ascii description */
	fread(&nchar, sizeof(short), (int) 1, fid);
	nchar = orderbyte(nchar);
	fseek(fid, nchar, SEEK_CUR);

	/* read rest of header */
	readshort(&(rfinfo->ncoils),  1, fid) ;
	readshort(&(rfinfo->res),     1, fid) ;
	readshort(&(rfinfo->npulses), 1, fid) ;
	fscanf(fid, "b1max:  %f\n", &(rfinfo->b1max));
	fscanf(fid, "gmax:   %f\n", &(rfinfo->gmax));

	fprintf(stderr, "\trfinfo.ncoils = %d \n", rfinfo->ncoils);
	fprintf(stderr, "\trfinfo.res    = %d \n", rfinfo->res);
	fprintf(stderr, "\trfinfo.npulses = %d \n", rfinfo->npulses);
/*
	fprintf(stderr, "\trfinfo.b1max  = %f \n", rfinfo->b1max);
	fprintf(stderr, "\trfinfo.gmax   = %f \n", rfinfo->gmax);
*/

	readshort(&(rfinfo->nparamsint16), 1, fid);
	readshort(rfinfo->paramsint16, rfinfo->nparamsint16, fid);

	rfinfo->npre  = rfinfo->paramsint16[0];
	rfinfo->rfres = rfinfo->paramsint16[1];

/*
	fprintf(stderr,"\tjfn_rf_readheader: npre, rfres = %d, %d\n", rfinfo->npre, rfinfo->rfres);

	rfinfo->npre  = 0;
	rfinfo->rfres = rfinfo->res;
*/

	readshort(&(rfinfo->nparamsfloat), 1, fid);
	for (i = 0; i < rfinfo->nparamsfloat; i++) {
		fscanf(fid, "%f\n", &(rfinfo->paramsfloat[i]));
		/* fprintf(stderr,"\tjfn_rf_readheader: rfinfo->paramsfloat[%d] = %f\n", i, rfinfo->paramsfloat[i]); */
	}

	rfinfo->dataoffset = ftell(fid);
  
	fclose (fid);
	return JFN_SUCCESS;
}


int readshort(short* i, int n, FILE *fid) {

	int j;

	fread(i, sizeof(short), n, fid);

	for (j=0; j<n; j++)
		i[j] = orderbyte(i[j]);

	return JFN_SUCCESS;
}


int readushort(unsigned short *i, int n, FILE *fid) {

	int j;

	fread(i, sizeof(unsigned short), n, fid);

	for (j=0; j<n; j++)
		i[j] = orderbyteu(i[j]);

	return JFN_SUCCESS;
}


int jfn_rf_readwaveforms(rfstruct *rfinfo, int ishard)
{
	FILE *fid;
	int p,c,i;   /* p = pulse number (0 or 1), c = coil #, i = waveform point index */
	int maxiamp = 32766;

	if ((fid = fopen (rfinfo->fname, "r")) == NULL)
		return(JFN_FAILURE);

	/* go past ascii and binary header */
	fseek(fid, rfinfo->dataoffset, SEEK_SET);

	/* fprintf(stderr,"jfn_rf_readwaveforms: ncoils, res, rfres, npre, dataoffset = %d, %d, %d, %d, %d \n", rfinfo->ncoils, rfinfo->res, rfinfo->rfres, rfinfo->npre, rfinfo->dataoffset); */

	/* read waveforms. Note that we're only using rf channel 1 in the psd itself. */
	for (p = 0; p < rfinfo->npulses; p++)  {
		for (c = 0; c < rfinfo->ncoils; c++) 
			readshort(rfinfo->rho[p][c],   (int) rfinfo->res, fid) ;
		for (c = 0; c < rfinfo->ncoils; c++) 
			readshort(rfinfo->theta[p][c], (int) rfinfo->res, fid) ;
		readshort(rfinfo->gx[p], (int) rfinfo->res, fid) ;
		readshort(rfinfo->gy[p], (int) rfinfo->res, fid) ;
		readshort(rfinfo->gz[p], (int) rfinfo->res, fid) ;
	}

	/* set to hard pulse if selected */
	if (ishard) {
		for (p = 0; p < rfinfo->npulses; p++)  {
			for (c = 0; c < rfinfo->ncoils; c++) {
				for (i = 0; i < rfinfo->res; i++) {
					rfinfo->rho[p][c][i]   = (short) maxiamp;
					rfinfo->theta[p][c][i] = (short) 0;
				}
			}
		}
	}

	/* make sure EOS bit of last point is set */
	for (p = 0; p < rfinfo->npulses; p++)  {
		for (c = 0; c < rfinfo->ncoils; c++) {
			if (!(rfinfo->rho[p][c][rfinfo->res-1] % 2))
				rfinfo->rho[p][c][rfinfo->res-1]++;
			if (!(rfinfo->theta[p][c][rfinfo->res-1] % 2))
				rfinfo->theta[p][c][rfinfo->res-1]++;
		}
		if (!(rfinfo->gx[p][rfinfo->res-1] % 2))
			rfinfo->gx[p][rfinfo->res-1]++;
		if (!(rfinfo->gy[p][rfinfo->res-1] % 2))
			rfinfo->gy[p][rfinfo->res-1]++;
		if (!(rfinfo->gz[p][rfinfo->res-1] % 2))
			rfinfo->gz[p][rfinfo->res-1]++;
	}

/*
#ifdef SIM
	for  (i = 0; i < rfinfo->res; i+=8) 
		fprintf(stderr,"\tjfn_rf_readwaveforms: rfinfo->gx[0][%d] = %d\n", i, rfinfo->gx[0][i]);
#endif
*/
  
	fclose (fid);
	return JFN_SUCCESS;
}


/* allocate memory for waveforms */
int jfn_rf_allocatemem(rfstruct *rfinfo) {

	int p,c;
	int npulses = rfinfo->npulses;

	rfinfo->rho   = (short***)AllocNode(sizeof(short**)*npulses);
	rfinfo->theta = (short***)AllocNode(sizeof(short**)*npulses);
	rfinfo->gx    = (short**)AllocNode(sizeof(short*)*npulses);
	rfinfo->gy    = (short**)AllocNode(sizeof(short*)*npulses);
	rfinfo->gz    = (short**)AllocNode(sizeof(short*)*npulses);

	for (p = 0; p < npulses; p++)  {
		rfinfo->rho[p]   = (short**)AllocNode(sizeof(short*)*rfinfo->ncoils);
		rfinfo->theta[p] = (short**)AllocNode(sizeof(short*)*rfinfo->ncoils);
		for (c=0; c < rfinfo->ncoils; c++) {
			rfinfo->rho[p][c]   = (short*)AllocNode(sizeof(short)*rfinfo->res);
			rfinfo->theta[p][c] = (short*)AllocNode(sizeof(short)*rfinfo->res);
		}
		rfinfo->gx[p] = (short*)AllocNode(sizeof(short)*rfinfo->res);
		rfinfo->gy[p] = (short*)AllocNode(sizeof(short)*rfinfo->res);
		rfinfo->gz[p] = (short*)AllocNode(sizeof(short)*rfinfo->res);
	}

	return JFN_SUCCESS;
}


/* free memory */
int jfn_rf_freemem(rfstruct *rfinfo) {

	int p,c;
	int npulses = rfinfo->npulses;

	for (p = 0; p < npulses; p++)  {
		FreeNode(rfinfo->gz[p]);
		FreeNode(rfinfo->gy[p]);
		FreeNode(rfinfo->gx[p]);
		for (c=0; c < rfinfo->ncoils; c++) {
			FreeNode(rfinfo->rho[p][c]);
			FreeNode(rfinfo->theta[p][c]);
		}
		FreeNode(rfinfo->theta[p]);
		FreeNode(rfinfo->rho[p]);
	}

	FreeNode(rfinfo->gz);
	FreeNode(rfinfo->gy);
	FreeNode(rfinfo->gx);
	FreeNode(rfinfo->theta);
	FreeNode(rfinfo->rho);

	return JFN_SUCCESS;
}

/* EOF */
