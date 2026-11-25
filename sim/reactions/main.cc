#include "Analysis.hh"
#include "Config.hh"
#include "ExperimentAssumptions.hh"
#include <filesystem>
#include "Telescope.hh"
#include "monad/monad.hxx"
#include "util/CLI.h"
#include "util/JSONParser.h"
#include "G4Box.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4GenericIon.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ProcessManager.hh"
#include "G4RunManager.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "TRandom.h"
#include "TROOT.h"
#include <locale>
#include <set>

using namespace indicators;

namespace {
class EventCountPunctuation final : public std::numpunct<char> {
	char do_thousands_sep() const override { return '\''; }
	std::string do_grouping() const override { return "\3"; }
};

phy::Nucleus parse_ion(const std::string& name) {
	auto n = phy::Nucleus::get_ion(name, true, "reaction list");
	if(!n.valid() || n.A < n.Z || n.Q() != n.Z) ERROR("Expected fully stripped ion: {}", name);
	return n;
}

std::vector<Resonance> channels(const nlohmann::json& list, const std::vector<int>& ids,
	const std::vector<double>& probabilities) {
	if(ids.size() != probabilities.size()) ERROR("Supply one probability for each reaction ID");
	std::set<int> seen;
	double total = 0;
	std::vector<Resonance> result;
	for(size_t i = 0; i < ids.size(); ++i) {
		const auto key = std::to_string(ids[i]);
		if(ids[i] <= 0 || !list.contains(key) || !seen.insert(ids[i]).second)
			ERROR("Unknown or repeated reaction ID {}", ids[i]);
		const double p = probabilities[i];
		if(!mnd::is_non_negative(p) || (total += p) > 1.) ERROR("Invalid reaction probabilities (sum must be <= 1)");
		const auto& spec = list.at(key);
		Resonance r{ids[i], parse_ion(spec.at("mother")), {}, spec.at("excitation").get<double>() * keV,
			spec.at("width").get<double>() * keV, total};
		if(!mnd::is_non_negative(r.excitation, r.width)) ERROR("Invalid resonance energy or width in {}", key);
		unsigned a = 0, z = 0;
		for(const auto& name: spec.at("daughters")) {
			const auto ion = parse_ion(name);
			a += ion.A;
			z += ion.Z;
			r.daughters.push_back(ion);
		}
		if(a != r.mother.A || z != r.mother.Z || r.daughters.size() < 2 || r.daughters.size() > 18)
			ERROR("Invalid daughter list for {}", key);
		std::stable_sort(r.daughters.begin(), r.daughters.end(), [](const auto& lhs, const auto& rhs) {
			return lhs.A > rhs.A;
		});
		if(!std::all_of(r.daughters.begin()+1, r.daughters.end(), [](const auto& n) { return n.A == 1 && n.Z == 1; }))
			ERROR("Expected a fragment and emitted protons in {}", key);
		result.push_back(std::move(r));
	}
	return result;
}

class Geometry final : public G4VUserDetectorConstruction {
	const transport::Config& telescope;
	double width_x, width_y, x, y;
public:
	double thickness;
	Geometry(const transport::Config& t, const nlohmann::json& target) : telescope(t),
		width_x(target.at("width_x")), width_y(target.at("width_y")), x(target.at("dx")), y(target.at("dy")),
		thickness(target.at("thickness").get<double>() * g/cm2 / (mnd::assume::s2::target_density*g/cm3)) {
		if(!mnd::is_positive(width_x, width_y, thickness)) ERROR("Invalid Be target dimensions");
	}

	G4VPhysicalVolume* Construct() override {
		auto* nist = G4NistManager::Instance();
		auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
		auto* world = new G4LogicalVolume(new G4Box("World", cfg::worldHalfXY, cfg::worldHalfXY,
			std::max(telescope.planes.back().z_mm, thickness) + 100*mm), vacuum, "World");
		auto* pv = new G4PVPlacement(nullptr, {}, world, "World", nullptr, false, 0, true);
		auto* be = nist->BuildMaterialWithNewDensity("ReactionBe", "G4_Be", (mnd::assume::s2::target_density*g/cm3));
		auto* target = new G4LogicalVolume(new G4Box("Target", width_x/2, width_y/2, thickness/2), be, "Target");
		new G4PVPlacement(nullptr, {x,y,-thickness/2}, target, "Target", world, false, 0, true);
		transport::place_telescope(telescope, world, 48., 0., false);
		return pv;
	}
};

class Physics final : public G4VModularPhysicsList {
	ReactionEvent& event;
public:
	explicit Physics(ReactionEvent& e) : event(e) {
		RegisterPhysics(new G4EmStandardPhysics_option4(0));
		SetDefaultCutValue(cfg::defaultProductionCut);
	}
	void ConstructProcess() override {
		G4VModularPhysicsList::ConstructProcess();
		G4GenericIon::Definition()->GetProcessManager()->AddDiscreteProcess(new Decay(event));
	}
	void SetCuts() override {
		SetCutsWithDefault();
		SetCutValue(cfg::deltaElectronRangeCut, "e-");
	}
};

class Beam final : public G4VUserPrimaryGeneratorAction {
	G4ParticleGun gun{1};
	const phy::Nucleus primary;
	const std::vector<Resonance>& choices;
	ReactionEvent& event;
	double thickness;
public:
	Beam(phy::Nucleus n, double energy, const std::vector<Resonance>& c, ReactionEvent& e, double t)
		: primary(n), choices(c), event(e), thickness(t) {
		gun.SetParticlePosition({0,0,-t-cfg::beamStartGap});
		gun.SetParticleMomentumDirection({0,0,1});
		gun.SetParticleEnergy(energy*n.A);
	}
	void GeneratePrimaries(G4Event* ge) override {
		event = {};
		const double draw = G4UniformRand();
		for(const auto& r: choices) {
			if(draw >= r.cumulative_probability) continue;
			event.channel = &r;
			event.requested_z = -thickness * G4UniformRand();
			event.measurements.resize(r.daughters.size());
			for(auto& track: event.measurements) track.fill(NAN);
			break;
		}
		gun.SetParticleDefinition(particle(primary));
		gun.GeneratePrimaryVertex(ge);
	}
};

class Crossings final : public G4UserSteppingAction {
	const transport::Config& geometry;
	ReactionEvent& event;
public:
	Crossings(const transport::Config& g, ReactionEvent& e) : geometry(g), event(e) {}
	void UserSteppingAction(const G4Step* step) override {
		auto* track = step->GetTrack();
		auto* id = dynamic_cast<DaughterID*>(track->GetUserInformation());
		if(!id) return;
		auto* pre = step->GetPreStepPoint();
		if(pre->GetStepStatus() != fGeomBoundary || pre->GetPhysicalVolume()->GetName() != "SiliconDown"
			|| pre->GetMomentumDirection().z() <= 0) return;
		const int i = pre->GetPhysicalVolume()->GetCopyNo();
		auto& measurement = event.measurements[id->index][i];
		const auto& plane = geometry.planes[i];
		const auto& pos = pre->GetPosition();
		// Fit in aligned coordinates: the strip offset cancels when converting back to world position.
		if(!std::isfinite(measurement) && std::abs(pos.z() - plane.z_mm) < 1e-6*mm)
			measurement = plane.a*pos.x() + plane.b*pos.y();
	}
};

class EndEvent final : public G4UserEventAction {
	ReactionEvent& event;
	Analysis& analysis;
	const u64 total;
	const u64 update_every;
	
	ProgressBar bar {
		option::BarWidth{55},
		option::Start{"["},
		option::Fill{"="},
		option::Lead{">"},
		option::Remainder{" "},
		option::End{"]"},
		option::PostfixText{fmt::format(std::locale(std::locale::classic(), new EventCountPunctuation),
			"Simulating reactions ({:L})", total)},
		option::ForegroundColor{Color::yellow},
		option::ShowPercentage{true},
		option::ShowElapsedTime{true},
		option::ShowRemainingTime{true},
		option::FontStyles{std::vector{FontStyle::bold}}
	};

public:
	EndEvent(ReactionEvent& e, Analysis& a, u64 count)
		: event(e), analysis(a), total(count), update_every((count + 99) / 100) {}

	void BeginOfEventAction(const G4Event* current) override {
		if(current->GetEventID() == 0) bar.set_progress(0);
	}

	void EndOfEventAction(const G4Event* current) override {
		analysis.record(event);
		mnd::PrintProgress(bar, current->GetEventID(), total, update_every);
	}
};
}

int main(int argc, char* argv[]) {
	CLI::App app{"Prompt decays in the 9Be target and propage the outward daughters out of the target, "
		"simulate the reconstruction in the FOOT telescope (assuming straggling). "
		"Is a direct simulation 'counterpart' to the quick/hit/physics program."};

	std::vector<int> ids;
	std::vector<double> probabilities;
	std::string list_path = REACTION_LIST, setup;
	u32 events = 10000;
	double energy = NAN;

	add_logged_option(app, "-w,--which", ids, "Reaction IDs in the JSON list")
		->delimiter(',')
		->required();
	add_logged_option(app, "-p,--probability", probabilities,
		"Probabilities per incident beam particle")
		->delimiter(',')
		->required();
	add_logged_option(app, "-l,--list", list_path, "Reaction list")
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-f,--setup", setup, "FOOT geometry")
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-e,--energy-per-u", energy,
		"Primary energy at target entrance [MeV/u]."
		"Default taken from the JSON file " REACTION_LIST
		)->check(CLI::PositiveNumber);
	add_logged_option(app, "-n,--events", events, "Number of incident particles");

	HistConfig hist_config;
	hist_config.add_options(app);

	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	if(test) return 0;
	
	srand(time(NULL));
	long seed = static_cast<long>(rand());

	const auto list = ParseJSON(list_path);
	const auto primary = parse_ion(list.at("secondary"));
	if(std::isnan(energy))
		energy = list.at("energy").get<double>();
	if(!mnd::is_positive(energy))
		ERROR("Expected finite positive beam energy in MeV/u");

	if(setup.empty()) setup = (std::filesystem::path(list_path).parent_path()
		/ list.at("target_file").get<std::string>()).string();

	const auto choices = channels(list.at("reaction"), ids, probabilities);
	const auto telescope = transport::load_config(setup, primary.to_string());

	ReactionEvent event {};
	gROOT->SetBatch(true);
	G4Random::setTheSeed(seed);
	gRandom->SetSeed(seed);
	
	const double nominal_reaction_energy =
		(energy + std::max(0.0, energy - mnd::assume::s2::loss_in_target)) / 2;
	Analysis analysis(telescope, hist_config, choices,
		phy::BetaGamma(primary, phy::EKin_t{nominal_reaction_energy}));
	auto manager = std::make_unique<G4RunManager>();
	auto* geometry = new Geometry(telescope, ParseJSON(setup).at("box").at("target"));
	
	show_console_cursor(false);
	manager->SetUserInitialization(geometry);
	manager->SetUserInitialization(new Physics(event));
	manager->SetUserAction(new Beam(primary, energy*MeV, choices, event, geometry->thickness));
	manager->SetUserAction(new Crossings(telescope, event));
	manager->SetUserAction(new EndEvent(event, analysis, events));
	manager->Initialize();
	manager->BeamOn(events);
	
	show_console_cursor(true);
	analysis.write();
}
