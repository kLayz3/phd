#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;

#pragma link C++ enum Orientation+;

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
#pragma link C++ class SCIQDCPedestal+;
#pragma link C++ class SCIMeanQDC+;
#pragma link C++ class SCIDEIntoQConverter+;
#pragma link C++ class SCIParam+;

#pragma link C++ class FRSIdParam+;
#pragma link C++ class FRSTargetParam+;
#pragma link C++ class FRSToFSingle+;
#pragma link C++ class FRSToFParam+;
#pragma link C++ class RNFRSHit+;
#pragma link C++ class RNFRSHit::Id+;

#pragma link C++ class RNFOOTMap+;
#pragma link C++ class FOOTClusterFit+;
#pragma link C++ class RNFOOTCluster+;
#pragma link C++ class RNFOOTCal+;
#pragma link C++ class FMultiPoly+;
#pragma link C++ class FOOTAsicGainParam+;
#pragma link C++ class FOOTReferentADCMeasurement+;
#pragma link C++ class FOOTGainParam+;
#pragma link C++ class FOOTDeltaFFT+;
#pragma link C++ class FOOTDeltaParam+;
#pragma link C++ class FOOTParam+;

#pragma link C++ class FOOTQ+;
#pragma link C++ class FOOTHit+;
#pragma link C++ class ExpertTarget+;
#pragma link C++ class FOOTBoxParam+;
#pragma link C++ class RNFOOTPair+;
#pragma link C++ class RNFOOTTrack+;
#pragma link C++ class RNFOOTHit+;
#pragma link C++ class RNFOOTHit::Vertex+;

#endif
