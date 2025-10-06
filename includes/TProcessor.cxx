#include "TProcessor.h"
#include "libs.hh"

Int_t TProcessor::Write() { 
	WARN("Write called but from the base class (no-op: Nothing is written). \n\tForgot to override it in derived class?\n");
	return 0; 
}
