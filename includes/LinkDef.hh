#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ class RNSciMap+;
#pragma link C++ class RNSciMap::Measurement+;
#pragma link C++ class RNTPCMap+;
#pragma link C++ class RNTPCMap::Measurement+;
#pragma link C++ class RNMUSICMap<8>+;
#pragma link C++ class RNFRSMap+;

#pragma link C++ class Scaler<32>+;
#pragma link C++ class RNTrigMap+;

#pragma link C++ class RNSciCal+;
#pragma link C++ class RNSciCal::Measurement+;
#pragma link C++ class RNTPCCal+;
#pragma link C++ class RNTPCCal::Measurement+;
#pragma link C++ class RNFRSCal+;
#pragma link C++ class TPCParam+;
#pragma link C++ class SCIParam+;

#pragma link C++ class RNFOOTMap++;
#pragma link C++ class RNFOOTCluster+;
#pragma link C++ class RNFOOTCal+;
#pragma link C++ class FOOTParam+;

#endif
