#include <iostream>
#include <cstdio>      // C-style file i/o
#include <string>
#include <stdexcept>
#include "toppe.h"
using namespace std;

// Class constructor
Module::Module(const char* s, int t) : type(t) {
	//type = t;

	try {
		readHeader(s);
		allocateMemory();

		try {
			readWaveforms(s);
		} catch (...) {
			cout << "Failed to read waveforms from " << s << endl;
			status = false;
		}

		status = true;
	} catch (std::exception& e){
		cout << "Exception: " << e.what() << endl;
	}
}

// Class destructor
Module::~Module() {
	if (status) {
		clearMemory();
	}
}
  
// Read header.
// Use C-style I/O so code can be copied easily from EPIC
void Module::readHeader(const char* fname) {

	FILE *fid;
	int p,c,i;   // p = pulse number (0 or 1), c = coil #, i = waveform point index

	if ((fid = fopen(fname, "r")) == NULL) {
		throw invalid_argument("Failed to open " + string(fname));
	}

	// go past ascii description 
	Module::readShort(&nascii, 1, fid);
	fseek(fid, nascii, SEEK_CUR);

	// read binary header
	readShort(&ncoils, 1, fid) ;
	readShort(&npts, 1, fid) ;
	readShort(&npulses, 1, fid) ;
	fscanf(fid, "b1max:  %f\n", &b1max);
	fscanf(fid, "gmax:  %f\n", &gmax);

	Module::readShort(&nparamsint16, 1, fid) ;
	Module::readShort(paramsint16, nparamsint16, fid);

	Module::readShort(&nparamsfloat, 1, fid) ;
	for (i = 0; i < nparamsfloat; i++) {
		fscanf(fid, "%f\n", &paramsfloat[i]);
	}

	dataoffset = ftell(fid);

	if (debug) {
		cout << "nascii: " << nascii << endl;
		cout << "ncoils: " << ncoils << endl;
		cout << "npts: " << npts << endl;
		cout << "npulses: " << npulses << endl;
		cout << "b1max: " << b1max << endl;
		cout << "gmax: " << gmax << endl;
		cout << "nparamsfloat: " << nparamsfloat << endl;
		cout << "nparamsint16: " << nparamsint16 << endl;
		cout << "dataoffset: " << dataoffset << endl;
	}

	fclose(fid);
}

void Module::readWaveforms(const char* fname) {

	int p, c, i;
	FILE *fid;

	fid = fopen(fname, "r");

	// go past header
	fseek(fid, dataoffset, SEEK_CUR);

	// read waveforms
	for (p = 0; p < npulses; p++)  {
		for (c = 0; c < ncoils; c++) {
			readShort(rho[p][c], (int) npts, fid) ;
		}
		for (c = 0; c < ncoils; c++) {
			readShort(theta[p][c], (int) npts, fid) ;
		}
		readShort(gx[p], (int) npts, fid) ;
		readShort(gy[p], (int) npts, fid) ;
		readShort(gz[p], (int) npts, fid) ;
	}

	fclose(fid);

	if (debug) {
		for  (i = npts-100; i < npts; i+=16) {
			cout << "gz[" << i << "]: " << gz[0][i]/(1.0*maxiamp)*gmax << endl;
		}
	}
}

void Module::allocateMemory() {
	int p, c; // counters for pulses, coils

	rho   = new short**[npulses];
	theta = new short**[npulses];
	gx    = new short*[npulses];
	gy    = new short*[npulses];
	gz    = new short*[npulses];

	for (p = 0; p < npulses; p++)  {
		rho[p]   = new short*[ncoils];
		theta[p] = new short*[ncoils];
		for (c=0; c < ncoils; c++) {
			rho[p][c]   = new short[npts];
			theta[p][c] = new short[npts];
		}
		gx[p] = new short[npts];
		gy[p] = new short[npts];
		gz[p] = new short[npts];
	}
}

void Module::clearMemory() {
	int p, c; // counters for pulses, coils

	for (p = 0; p < npulses; p++)  {
		for (c=0; c < ncoils; c++) {
			delete[] rho[p][c];
			delete[] theta[p][c];
		}
		delete[] rho[p];
		delete[] theta[p];
		delete[] gx[p];
		delete[] gy[p];
		delete[] gz[p];
	}

	delete[] rho;
	delete[] theta;
	delete[] gx;
	delete[] gy;
	delete[] gz;
}


// Read short int array of length n
void Module::readShort(short* sa, int n, FILE *fid) {
	int j;
	fread(sa, sizeof(short), n, fid);
	for (j=0; j<n; j++) {
		sa[j] = Module::swapByte(sa[j]);  
	}
}

// convert from Big to Little endian
short Module::swapByte(short in)
{
  char* sw;
  char tmp;
  sw = (char*)&in;
  tmp=sw[0];
  sw[0] = sw[1];
  sw[1] = tmp;
  return in;
}
