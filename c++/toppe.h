#ifndef TOPPE_H
#define TOPPE_H

#include <fstream>
#include <cstdio>
#include <string>
using namespace std;

#define FAILURE 0
#define SUCCESS 1

class Module {
private:
	bool debug = true;        // If true, prints some values from the .mod file to screen
	bool status = false;      // set to true once this Module has successfully loaded a .mod file

	short nascii;             // number of characters in ascii header
	short dataoffset;         // total header size (including ASCII and binary parts)
	short nparamsint16;       // # of int16 header parameters
	short nparamsfloat;       // # of float header parameters

	void  readShort(short*, int, FILE* );
	short swapByte(short);

public:
	short    ncoils;             // number of coils/channels. Useful in the future for, e.g., parallel transmit systems.
	short    npts;               // number of points in waveform
	short    npulses;            // number of different waveforms in .wav file 
	float    b1max;              // RF amplitude waveform scaling (Gauss). Used to convert to physical units.
	float    gmax;               // Gradient waveform scaling (Gauss/cm). Used to convert physical units.
	short    paramsint16[32];    // int16 parameters 
	float    paramsfloat[32];    // float parameters

	int maxiamp = 32766;        /* Max waveform amplitude; In TOPPE, waveforms are saved in int16 format */

	/* 
	RF and gradient waveform arrays. To convert these to physical units:
	rho   = rho/maxiamp*b1max   (Gauss)
	theta = theta/maxiamp*pi    (Radians)
	gx    = gx/maxiamp*gmax     (Gauss/cm)
	gy    = gy/maxiamp*gmax     (Gauss/cm)
	gz    = gz/maxiamp*gmax     (Gauss/cm)
	*/
	short*** rho;               // RF magnitude waveforms (one or more, possibly multi-coil). rho[npulses][ncoils][npts]
	short*** theta;             // RF phase waveform(s).  theta[npulses][ncoils][npts]
	short**  gx;                // x gradient waveform(s). gx[npulses][npts]
	short**  gy;
	short**  gz;

	// Information from modules.txt
	int type;                  // 0: gradients only; 1: RF excitation module; 2: DAQ module
	int duration;              // If duration > minimum module duration, a delay is added to end of module (allows variable TR)

	// Constructor and destructor
	Module(const char* s, int type);
	~Module();                 // for clearing memory

	void readHeader(const char*);
	void readWaveforms(const char*);
	void allocateMemory();
	void clearMemory();
};

#endif
