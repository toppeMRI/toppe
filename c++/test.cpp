#include <iostream>
#include <string>
#include "toppe.h"
using namespace std;

int main() {
	int type;                  // 0: gradients only; 1: RF excitation module; 2: DAQ module
	type = 1;
	Module ex = Module("tipdown.mod", type);
	Module no = Module("existsnot.mod", type);
}
