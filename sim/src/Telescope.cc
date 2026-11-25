#include "Telescope.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

void transport::place_telescope(const Config& config, G4LogicalVolume* world_lv,
	double half_width, double origin_z, bool use_vacuum) {
	auto* nist = G4NistManager::Instance();
	auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
	double before_si = 0;
	for(const auto& layer: config.materials) {
		if(layer.name == "G4_Si") {
			before_si += layer.thickness_mm / 2;
			break;
		}
		before_si += layer.thickness_mm;
	}
	for(size_t i=0; i<N; ++i) {
		const auto& plane = config.planes[i];
		double cursor = plane.z_mm - before_si;
		for(const auto& layer: config.materials) {
			auto* material = use_vacuum ? vacuum : nist->FindOrBuildMaterial(layer.name);
			if(!material) throw std::runtime_error("Unknown Geant4 material: " + layer.name);
			const bool silicon = layer.name == "G4_Si";
			const int pieces = silicon ? 2 : 1;
			for(int part=0; part<pieces; ++part) {
				const double t = layer.thickness_mm/pieces;
				const std::string name = silicon ? (part ? "SiliconDown" : "SiliconUp") : "Passive";
				auto* lv = new G4LogicalVolume(new G4Box(name,half_width*mm,half_width*mm,t*mm/2),material,name);
				new G4PVPlacement(nullptr,
					{plane.a*plane.offset_mm*mm,plane.b*plane.offset_mm*mm,(cursor+t/2-origin_z)*mm},
					lv,name,world_lv,false,static_cast<int>(i),true);
				cursor += t;
			}
		}
	}
}
