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

  $Id: cores.h,v 1.13 2017/08/28 22:15:55 jfnielse Exp $
 */

#ifndef CORES_H
#define CORES_H

typedef struct {
	int ncores;           /* number of modules (cores) */
	char** wavfiles; 
	int* dur;     /* msec. 0 corresponds to minimum realizable duration. */
	int* hasRF;
	int* hasDAQ;
} corestruct;

int cores_getinfo(char *fname, corestruct* coredefinfo);

int cores_getloophdr(char *fname, int* loophdr);

int cores_readloop(char *fname, int* looparr);

#endif
