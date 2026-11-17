#include "Config.hh"
#include "util/MPhysics.h"
#include "util/CLI.h"
#include "monad/accessories/typedefs.hxx"

#include "G4AnalysisManager.hh"
#include "G4Box.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace transport {
struct Options {
	std::string setup = FOOT_PARAM_DIR "/foot_setup.json";
	std::string cov = FOOT_PARAM_DIR "/foot_cov.json";
	std::string ion = "9C", output = "foot_transport.root";
	double energy = 0, half_width = 48, x = 0, y = 0, ax = 0, ay = 0;
	double beam_sigma_mm = 0, angle_sigma_mrad = 0, max_step_mm = 0.015;
	int events = 1000, threads = 1;
	long seed = 12345;
	bool check = false, vacuum = false;
	int A = 9, Z = 6;
};

class Geometry final : public G4VUserDetectorConstruction {
	const Config& config;
	const Options& options;
public:
	double centre, entrance;
	Geometry(const Config& c, const Options& o) : config(c), options(o) {
		centre = (c.planes.front().z_mm + c.planes.back().z_mm)/2;
		double before_si = 0;
		for(const auto& layer: c.materials) {
			if(layer.name == "G4_Si") {before_si += layer.thickness_mm/2; break;}
			before_si += layer.thickness_mm;
		}
		entrance = c.planes.front().z_mm - before_si - 1;
	}
	G4VPhysicalVolume* Construct() override {
		auto* nist = G4NistManager::Instance();
		auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
		double thickness = 0, before_si = 0;
		bool found_si = false;
		for(const auto& layer: config.materials) {
			thickness += layer.thickness_mm;
			if(!found_si) {
				before_si += layer.thickness_mm * (layer.name == "G4_Si" ? 0.5 : 1);
				found_si = layer.name == "G4_Si";
			}
		}
		for(size_t i=1; i<N; ++i)
			if(config.planes[i].z_mm-config.planes[i-1].z_mm <= thickness)
				throw std::runtime_error("Detector material stacks overlap");
		const double world_xy = (options.half_width + std::max(std::abs(options.x), std::abs(options.y)) + 100)*mm;
		const double world_z = ((config.planes.back().z_mm-config.planes.front().z_mm)/2 + thickness + 100)*mm;
		auto* world_lv = new G4LogicalVolume(new G4Box("World",world_xy,world_xy,world_z),vacuum,"World");
		auto* world = new G4PVPlacement(nullptr,{},world_lv,"World",nullptr,false,0,true);
		for(size_t i=0; i<N; ++i) {
			const auto& plane = config.planes[i];
			double cursor = plane.z_mm - before_si;
			for(const auto& layer: config.materials) {
				auto* material = options.vacuum ? vacuum : nist->FindOrBuildMaterial(layer.name);
				if(!material) throw std::runtime_error("Unknown Geant4 material: " + layer.name);
				const bool silicon = layer.name == "G4_Si";
				const int pieces = silicon ? 2 : 1;
				for(int part=0; part<pieces; ++part) {
					const double t = layer.thickness_mm/pieces;
					const std::string name = silicon ? (part ? "SiliconDown" : "SiliconUp") : "Passive";
					auto* lv = new G4LogicalVolume(new G4Box(name,options.half_width*mm,options.half_width*mm,t*mm/2),material,name);
					new G4PVPlacement(nullptr,
						{plane.a*plane.offset_mm*mm,plane.b*plane.offset_mm*mm,(cursor+t/2-centre)*mm},
						lv,name,world_lv,false,static_cast<int>(i),true);
					cursor += t;
				}
			}
		}
		return world;
	}
};

class Physics final : public G4VModularPhysicsList {
public:
	Physics() {
		SetVerboseLevel(0);
		SetDefaultCutValue(0.1*mm);
		RegisterPhysics(new G4EmStandardPhysics_option4(0));
		// Step limitation is included explicitly below via G4StepLimiterPhysics.
	}
	void SetCuts() override {
		SetCutsWithDefault();
		SetCutValue(10*mm, "e-");
	}
};

struct Crossing {
	bool present = false;
	double x, y, z, ax, ay, energy, measurement;
};
struct EventData {
	std::array<Crossing,N> planes{};
	double x, y, ax, ay;
};

class Generator final : public G4VUserPrimaryGeneratorAction {
	const Options& o;
	const Geometry& geometry;
	EventData& data;
	G4ParticleGun gun{1};

public:
	Generator(const Options& opts,const Geometry& g,EventData& d) : o(opts),geometry(g),data(d) {}
	void GeneratePrimaries(G4Event* event) override {
		data.x = o.x + G4RandGauss::shoot(0,o.beam_sigma_mm);
		data.y = o.y + G4RandGauss::shoot(0,o.beam_sigma_mm);
		data.ax = o.ax + G4RandGauss::shoot(0,o.angle_sigma_mrad*0.001);
		data.ay = o.ay + G4RandGauss::shoot(0,o.angle_sigma_mrad*0.001);
		auto* particle = o.A == 1 && o.Z == 1 ? G4Proton::ProtonDefinition()
			: G4IonTable::GetIonTable()->GetIon(o.Z,o.A,0);
		gun.SetParticleDefinition(particle);
		gun.SetParticleCharge(o.Z*eplus);
		gun.SetParticleEnergy(o.A*o.energy*MeV);
		gun.SetParticlePosition({data.x*mm,data.y*mm,(geometry.entrance-geometry.centre)*mm});
		gun.SetParticleMomentumDirection(G4ThreeVector(data.ax,data.ay,1).unit());
		gun.GeneratePrimaryVertex(event);
	}
};

class EventAction final : public G4UserEventAction {
public:
	EventData data;
	void BeginOfEventAction(const G4Event*) override {data.planes = {};}
	void EndOfEventAction(const G4Event* event) override {
		auto* a = G4AnalysisManager::Instance();
		a->FillNtupleIColumn(0,event->GetEventID());
		a->FillNtupleDColumn(1,data.x); a->FillNtupleDColumn(2,data.y);
		a->FillNtupleDColumn(3,data.ax); a->FillNtupleDColumn(4,data.ay);
		int column = 5;
		for(const auto& p: data.planes) {
			a->FillNtupleIColumn(column++,p.present ? 1 : 0);
			const double nan = std::numeric_limits<double>::quiet_NaN();
			for(double value: {p.x,p.y,p.z,p.ax,p.ay,p.energy,p.measurement})
				a->FillNtupleDColumn(column++,p.present ? value : nan);
		}
 		a->AddNtupleRow();
	}
};

class Stepping final : public G4UserSteppingAction {
	EventData& data;
	const Config& config;
	const Geometry& geometry;
public:
	Stepping(EventData& d,const Config& c,const Geometry& g) : data(d),config(c),geometry(g) {}
	void UserSteppingAction(const G4Step* step) override {
		const auto* track = step->GetTrack();
		if(track->GetParentID() != 0 || track->GetTrackID() != 1) return;
		const auto* post = step->GetPostStepPoint();
		const auto* pre_volume = step->GetPreStepPoint()->GetPhysicalVolume();
		const auto* post_volume = post->GetPhysicalVolume();
		if(post->GetStepStatus() != fGeomBoundary || !pre_volume || !post_volume) return;
		// Splitting silicon at its midpoint makes the measurement a real transport boundary.
		if(pre_volume->GetName() != "SiliconUp" || post_volume->GetName() != "SiliconDown") return;
		const int index = post_volume->GetCopyNo();
		if(index < 0 || index >= static_cast<int>(N) || pre_volume->GetCopyNo() != index) return;
		auto& p = data.planes[index];
		if(p.present) return;
		const auto direction = post->GetMomentumDirection();
		if(direction.z() <= 0) return;
		const auto position = post->GetPosition()/mm;
		const auto& plane = config.planes[index];
		p = {true,position.x(),position.y(),position.z()+geometry.centre,
			direction.x()/direction.z(),direction.y()/direction.z(),post->GetKineticEnergy()/MeV,
			plane.a*position.x()+plane.b*position.y()-plane.offset_mm};
	}
};

class RunAction final : public G4UserRunAction {
	const Options& options;
public:
	explicit RunAction(const Options& o) : options(o) {
		auto* a = G4AnalysisManager::Instance();
		a->SetVerboseLevel(0);
		a->SetDefaultFileType("root");
		a->SetNtupleMerging(o.threads > 1);
		a->CreateNtuple("crossings","Unsmeared primary midplane crossings; detector index order");
		a->CreateNtupleIColumn("event");
		for(const char* name: {"beam_x_mm","beam_y_mm","beam_ax","beam_ay"}) a->CreateNtupleDColumn(name);
		for(size_t i=0; i<N; ++i) {
			const auto prefix = "d" + std::to_string(i) + "_";
			a->CreateNtupleIColumn(prefix+"present");
			for(const char* name: {"x_mm","y_mm","z_mm","ax","ay","kinetic_MeV","measurement_mm"})
				a->CreateNtupleDColumn(prefix+name);
		}
		a->FinishNtuple();
	}
	void BeginOfRunAction(const G4Run*) override {
		if(!G4AnalysisManager::Instance()->OpenFile(options.output))
			throw std::runtime_error("Cannot open transport output");
	}
	void EndOfRunAction(const G4Run*) override {
		auto* a = G4AnalysisManager::Instance(); a->Write(); a->CloseFile();
	}
};

class Actions final : public G4VUserActionInitialization {
	const Config& config;
	const Options& options;
	const Geometry& geometry;
public:
	Actions(const Config& c,const Options& o,const Geometry& g) : config(c),options(o),geometry(g) {}
	void BuildForMaster() const override {SetUserAction(new RunAction(options));}
	void Build() const override {
		auto* event = new EventAction;
		SetUserAction(event);
		SetUserAction(new Generator(options,geometry,event->data));
		SetUserAction(new RunAction(options));
		SetUserAction(new Stepping(event->data,config,geometry));
	}
};
}

#include "G4StepLimiterPhysics.hh"
#include "G4UserLimits.hh"

int main(int argc,char** argv) {
	using namespace transport;
	Options o;
	CLI::App app{"FOOT electromagnetic telescope transport; no reaction or detector smearing"};
	add_logged_option(app,"--setup", o.setup, "Detector geometry and material parameter file")
		->check(CLI::ExistingFile);
	add_logged_option(app,"--cov", o.cov, "Measured charge/position residual parameter file")
		->check(CLI::ExistingFile);
	add_logged_option(app,"--ion", o.ion, "Isotope parsed by phy::Nucleus, e.g. 1H, 4He, 9C");
	add_logged_option(app,"--energy-per-u", o.energy, "Kinetic energy at telescope entrance [MeV/u]")
		->required()
		->check(CLI::PositiveNumber);
	add_logged_option(app,"--output", o.output, "ROOT output; companion .json records model inputs");
	add_logged_option(app,"--events", o.events, "Number of primary events")
		->check(CLI::PositiveNumber);
	add_logged_option(app,"--threads", o.threads, "Number of Geant4 workers (requires an MT build)")
		->check(CLI::PositiveNumber);
	add_logged_option(app,"--seed", o.seed, "Random-number seed")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--half-width-mm", o.half_width,
		"Square detector half-width (geometry approximation)")
		->check(CLI::PositiveNumber);
	add_logged_option(app,"--x-mm", o.x, "Beam entrance x [mm]");
	add_logged_option(app,"--y-mm", o.y, "Beam entrance y [mm]");
	add_logged_option(app,"--ax", o.ax, "Beam entrance dx/dz");
	add_logged_option(app,"--ay",o.ay, "Beam entrance dy/dz");
	add_logged_option(app,"--beam-sigma-mm", o.beam_sigma_mm, "Gaussian beam position RMS in each axis [mm]")
		->check(CLI::NonNegativeNumber);
	add_logged_option(app,"--angle-sigma-mrad", o.angle_sigma_mrad,
		"Gaussian beam slope RMS in each axis [mrad]")
		->check(CLI::NonNegativeNumber);
	add_logged_option(app, "--max-step-mm", o.max_step_mm,
		"Maximum silicon/material step; vary for sanity checks")
		->check(CLI::PositiveNumber);
	add_logged_flag(app,"--vacuum", o.vacuum, "Replace material with vacuum for exact straight-line validation");
	add_logged_flag(app,"--check-config", o.check,
		"Parse and print configuration without running or writing output");

	CLI11_PARSE(app,argc,argv);
	
	auto nucleus = phy::Nucleus::get_ion(o.ion, true, "foot_transport --ion");
	if(!nucleus.valid() || nucleus.A < nucleus.Z || nucleus.Q() != nucleus.Z)
		ERROR("Transport requires a valid fully stripped ion.");
	o.A = static_cast<int>(nucleus.A);
	o.Z = static_cast<int>(nucleus.Z);
	nucleus.N_electrons.reset();
	o.ion = nucleus.to_string();

	for(double v: {o.energy,o.half_width,o.x,o.y,o.ax,o.ay,o.beam_sigma_mm,o.angle_sigma_mrad,o.max_step_mm})
		if(!std::isfinite(v)) ERROR("Beam and geometry options must be finite");
	
	auto config = load_config(o.setup,o.cov,o.ion);

	config.provenance["nucleus"] = {
		{"A",nucleus.A}, {"Z",nucleus.Z}, {"phy_mass_MeV", nucleus.mass()}
	};
	config.provenance["transport"] = {
		{"energy_per_u_MeV", o.energy},
		{"events",           o.events},
		{"threads",          o.threads},
		{"seed",             o.seed},
		{"half_width_mm",    o.half_width},
		{"vacuum",           o.vacuum},
		{"max_step_mm",      o.max_step_mm},
		{"beam_x_mm",        o.x},
		{"beam_y_mm",        o.y},
		{"beam_ax",          o.ax},
		{"beam_ay",          o.ay},
		{"beam_sigma_mm",    o.beam_sigma_mm},
		{"angle_sigma_mrad", o.angle_sigma_mrad},
		{"physics",          "EM option4; no hadronic reactions or decay; 10 mm electron production range cut"},
		{"scope",            "telescope material only; no target, air, or supports unless specified in composition"}
	};

	if(o.check) {
		std::cout << config.provenance.dump(4) << '\n'; 
		return 0;
	}
	
	if(std::filesystem::path(o.output).extension() != ".root")
		ERROR("Output must have .root extension");
	if(std::filesystem::exists(o.output) || std::filesystem::exists(o.output+".json"))
		ERROR("Output already exists; choose a new --output path");

#ifndef G4MULTITHREADED
	if(o.threads != 1)
		throw std::runtime_error("This Geant4 installation is serial-only; use --threads 1 or an MT build");
#endif

	G4Random::setTheSeed(o.seed);
	auto manager = std::unique_ptr<G4RunManager>(
		G4RunManagerFactory::CreateRunManager(
			o.threads > 1 ? G4RunManagerType::MTOnly : G4RunManagerType::SerialOnly
		)
	);
	manager->SetNumberOfThreads(o.threads);

	auto* geometry = new Geometry(config,o);
	manager->SetUserInitialization(geometry);
	
	auto* physics = new Physics;

	physics->RegisterPhysics(new G4StepLimiterPhysics);
	manager->SetUserInitialization(physics);
	manager->SetUserInitialization(new Actions(config,o,*geometry));
	manager->Initialize();

	config.provenance["transport"]["actual_threads"] = manager->GetNumberOfThreads();

	auto* particle = o.A == 1 && o.Z == 1 ? G4Proton::ProtonDefinition()
		: G4IonTable::GetIonTable()->GetIon(o.Z,o.A,0);
	
	if(!particle)
		ERROR("Geant4 cannot construct the selected nucleus");
	
	config.provenance["nucleus"]["geant4_mass_MeV"] = particle->GetPDGMass()/MeV;
	
	// Logical volumes own neither the limits nor the configuration; retain limits until the run manager is destroyed.
	auto limits = std::make_unique<G4UserLimits>(o.max_step_mm*mm);
	for(auto* lv: *G4LogicalVolumeStore::GetInstance())
		if(lv->GetName() != "World") lv->SetUserLimits(limits.get());
	config.provenance["beam_entrance_z_mm"] = geometry->entrance;
	
	std::ofstream metadata(o.output + ".json");
	if(!metadata) ERROR("Cannot open transport metadata output");
	metadata << config.provenance.dump(4) << '\n'; metadata.close();
	
	if(!metadata) ERROR("Cannot write transport metadata");
	
	manager->BeamOn(o.events);
	manager.reset();
}
