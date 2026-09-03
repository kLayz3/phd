/* Main program: on the track multp. vs charge graph do different cuts and look
 * at the resulting spectra... */

#include "util/CLI.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TParameter.h"
#include "util/Geometry.h"
#include "util/MacroHelpers.h"
#include "util/PrettyHisto.h"

#include "TFOOTHitCont.h"
#include "TFRSHitCont.h"
#include "TFRSCalCont.h"
#include <cassert>
#include <cmath>

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

#define MND_PHYSICS_PARANOIA

enum class AtomicNumber { H, He, Li, Be, B, C };
enum class AngleType { all, p, he };

constexpr A2 elem_to_a2(AtomicNumber e) {
	switch(e) {
		case AtomicNumber::H:  return A2{0.6, 1.5};
		case AtomicNumber::He: return A2{1.5, 2.5};
		case AtomicNumber::Li: return A2{2.5, 3.5};
		case AtomicNumber::Be: return A2{3.5, 4.5};
		case AtomicNumber::B:  return A2{4.5, 5.4};
		case AtomicNumber::C:  return A2{5.4, 6.8};
	}
	return A2{};
}

struct SingleSelect {
	std::variant<AtomicNumber, A2> fragment;
	u32 index;
};
std::istream& operator>>(std::istream&, SingleSelect& );
std::ostream& operator<<(std::ostream&, const SingleSelect& );

using Select = std::vector<SingleSelect>;

int main(int argc, char* argv[]) {
	CLI::App app{"This program analyses the Hit level file and produces (hopefully) decent physics."};

	std::vector<std::string> fileName, info;
	double distance_cut = NAN;
	auto sci21_cut = mnd::make_filled_array<double,2>(NAN);
	auto sci22_cut = mnd::make_filled_array<double,2>(NAN);
	auto sci31_cut = mnd::make_filled_array<double,2>(NAN);
	size_t max_events = -1;

	Select selected{};
	auto angle_type = AngleType::all;
	std::vector<canvas::Extension> save = {};

	add_logged_option(app, "-f,--file", fileName, "Pass one or more file names, delimited by ','")
		->delimiter(',')
		->check(CLI::ReadPermissions);
	
	add_logged_option<DisplayDefault::No>(app, "-m,--max-events", max_events,
		"Max events (at most) taken from each ROOT file. Default: all entries.");
	add_logged_option(app, "-s,--select", selected, "Select which fragments to explicity gate upon, delimited by ':'")
		->delimiter(':')
		->type_name("ELEM,N | A,B,N");

	add_logged_option<DisplayDefault::No>(app, "-d,--distance", distance_cut,
		"Distance from vertex cut (all tracks must be below this threshold). Default no cut. "
		"NOTE: not assigning this field will default to simple angles between all FOOT tracks.")
		->check(CLI::Range(0.0, 2.0));
	add_logged_option<DisplayDefault::No>(app, "--sci21",sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22",sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31",sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');

	add_logged_option(app, "-t,--angle", angle_type, "Specify which type of formula to use for angular "
		"spectroscopic ρ coefficient of N-body decay:\n"
		"[he] : sqrt(Σ θ_i^2) angles between all heavy tracks and a single light one, (N-1) angles.\n"
		"[p]  : sqrt(Σ θ_i^2) angles between all light tracks and a single heavy one, (N-1) angles.\n" 
		"[all]: sqrt(Σ θ_i^2) combination of all the angles between all the tracks, N*(N-1)/2 angles.");
	add_logged_option(app, "-o,--save", save, "Save the resulting canvases as one or more extensions.")
		->delimiter(',');
	add_logged_option(app, "-i,--info", info,
		"Sequence of strings, such that the output will be saved descending these directories in autosave.")
		->delimiter(',');
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;
	if(fileName.size() == 0)
		ERROR("To continue, must supply at least one file name!\n");

	mnd::python::poke();

	TApplication rootApp("app", 0, 0);

	using namespace mnd::geom;

	double Cr, Cq, Ct, max_cost, max_cost_f;
	{
		const auto& fname = fileName.front();
		std::array<double, 3>* c;
		TParameter<double>* m;
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fname.c_str(), "READ");
		get_obj(f, c, "FOOT_cost_coeff");
		get_obj(f, m, "FOOT_max_cost");
		Cr = c->at(0); Cq = c->at(1); Ct = c->at(2);
		max_cost = m->GetVal();
		get_obj(f, m, "FOOT_max_cost_f");
		max_cost_f = m->GetVal();
	}

	/* Sanitize some CLI passed in arguments... */
	const u32 sum_n_tracks_required = std::accumulate(
		selected.begin(), selected.end(), 0, [](u32 sum, const auto& x) {
			return sum + x.index;
		}
	);
	if(sum_n_tracks_required == 1)
		ERROR("Cannot require just a single track to form a valid vertex. -s,--select option "
		      "must require at least two particles.\n");

	/* Convert the selected variant exclusively into the A2 = std::array<double, 2>. */
	for(auto& [key,value] : selected) {
		if(value == 0) 
			ERROR("Must require 1 or more tracks, cannot set to require 0 tracks. Stop that.");
		if(std::holds_alternative<AtomicNumber>(key)) {
			key = elem_to_a2( std::get<AtomicNumber>(key) );
		}
	}
	
	/* All cut entries now hold only A2 range. Sort them out, in descending charge.
	 * Convention from RNFOOTTrack::operator<(...) .
	 * Also use this time to do some sanity checks for the ranges. */
	std::sort(selected.begin(), selected.end(),
		/* Sort in descending order, according to the lower range boundary. */
		[](const auto& lhs, const auto& rhs) {
			if(std::get<1>(lhs.fragment)[0] > std::get<1>(lhs.fragment)[1]) {
				WARN("Range: "); std::cerr << std::get<1>(lhs.fragment) << " isn't valid (lhs >= rhs)?";
				throw std::runtime_error("Err");
			}
			if(std::get<1>(rhs.fragment)[0] > std::get<1>(rhs.fragment)[1]) {
				WARN("Range: "); std::cerr << std::get<1>(rhs.fragment) << " isn't valid (rhs >= rhs)?";
				throw std::runtime_error("Err");
			}
			return std::get<1>(lhs.fragment)[0] > std::get<1>(rhs.fragment)[0];
		}
	);
	if( std::adjacent_find(selected.begin(), selected.end(),
		[](const auto& lhs, const auto& rhs) {
			return std::get<1>(lhs.fragment)[0] < std::get<1>(rhs.fragment)[1];
			// If LHS elem's lower bound is smaller than RHS elem's upper bound, means their intervals overlap.
		}
	) != selected.end() ) {
		ERROR("Two intervals from the 'selected' (-s,--select) field are overlapping. Not allowed");
	}

	/* Sanitize the angle part. If selected isn't given by at least
	 * two windows, then angle cannot be clearly defined between light and heavy ones. */
	if(selected.size() < 1 and (angle_type == AngleType::p || angle_type == AngleType::he)) {
		angle_type = AngleType::all;
	}

	const bool vertex_dist_cut_given = std::isfinite(distance_cut);

	TH1P* h1_track_mult = new TH1P("Track multiplicity [unique tracks]@-1 means no FOOT in event.", kRed-1, 11, -1.5, 9.5);
	if(vertex_dist_cut_given and sum_n_tracks_required > 0)
		h1_track_mult->AppendToTitle(" 0 means event not viable.");

	TH2P* h2_q_vs_mult = new TH2P("Track charge [Q]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 40, 0,8);
	TH2P* h2_score_vs_mult = new TH2P("Track score [a.u.]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 300, 0, 50);
	TH2P* h2_track_distance = new TH2P("Track distances to vertex [mm]:Track I@smaller x-axis means larger Q particle",
		(vertex_dist_cut_given && sum_n_tracks_required>0)? sum_n_tracks_required + 1: 10,
		-0.5,
		(vertex_dist_cut_given && sum_n_tracks_required>0)? sum_n_tracks_required-0.5: 9.5,
		600, 0, vertex_dist_cut_given? (2*distance_cut): 30.0);

	TH1P* h1_track_angle = new TH1P("Track angles [mrad]@Between all tracks selected", kMagenta+1, 200, 0, 100);
	TH1P* h1_angle_ex = new TH1P (
		Form("#sqrt( #sum_{i=1}^{%s} #theta_{i|heavy}^{2} ) [mrad]@#rho angle",
			(sum_n_tracks_required>0 && vertex_dist_cut_given) ? std::to_string(sum_n_tracks_required-1).c_str(): "??"),
		kCyan-9, 300, 0, 300
	);
	if(sum_n_tracks_required == 2 and angle_type == AngleType::p)
		(*h1_angle_ex)->GetXaxis()->SetTitle("#theta(p,frag) [mrad]");

	TH2P* h2_vertex_z = new TH2P("#rho angle [mrad]:Vertex z [mm]@Traced by the FOOT", 160, -80, 80, 100,0,100);
	
	auto* h1_sci21 = new TH1P("SCI21 QDC mean [QDC units]", 0xCB00CB_c, 500, 300, 4000);
	auto* h1_sci22 = new TH1P("SCI22 QDC mean [QDC units]", 0x0070DD_c, 500, 300, 4000);
	auto* h1_sci31 = new TH1P("SCI31 QDC mean [QDC units]", 0x009B2F_c, 500, 300, 4000);
	auto* h1_sci21_cut  = new TH1P("((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, 500, 300, 4000);
	auto* h1_sci22_cut  = new TH1P("((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, 500, 300, 4000);
	auto* h1_sci31_cut  = new TH1P("((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0x7DE69D_c, 500, 300, 4000);
	auto* h1_sci21_cut2 = new TH1P("((h1_cut2)) SCI21 QDC mean [QDC units]@With cut and FOOT selection", 0x890389_c, 500, 300, 4000);
	auto* h1_sci22_cut2 = new TH1P("((h1_cut2)) SCI22 QDC mean [QDC units]@With cut and FOOT selection", 0x6180FD_c, 500, 300, 4000);
	auto* h1_sci31_cut2 = new TH1P("((h1_cut2)) SCI31 QDC mean [QDC units]@With cut and FOOT selection", 0x7DE69D_c, 500, 300, 4000);

	show_console_cursor(false);
	
	for(size_t i{0}; i < fileName.size(); ++i) {
		const auto& fname = fileName[i];
		auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTHit>("FOOT");
		auto frs = model->MakeField<RNFRSHit>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);
		const size_t nentries = ( (max_events < ntuple->GetNEntries()) ? max_events : ntuple->GetNEntries() );
		ProgressBar bar {
			option::BarWidth{50},
			option::Start{"["},
			option::Fill{"="},
			option::Lead{">"},
			option::Remainder{" "},
			option::End{"]"},
			option::PostfixText{mnd::msg("Analysis (per event: %s)", fname.c_str())},
			option::ForegroundColor{Color::yellow},
			option::ShowPercentage{true},
			option::ShowElapsedTime{true},
			option::ShowRemainingTime{true},
			option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
		};
		WARN("Proceeding with file [%zu/%zu]: \'%s\'. Entries: [%'zu]\n", i+1, fileName.size(), fname.c_str(), nentries);
		
		for(size_t entryId{0}; entryId < nentries; ++entryId ) {
			mnd::PrintProgress(bar, entryId, nentries, 500, mnd::dancer0, 0.33);

			ntuple->LoadEntry(entryId);

			const auto& sci21 = frs->cal.sci[0];
			const auto& sci22 = frs->cal.sci[1];
			const auto& sci31 = frs->cal.sci[2];
			
			if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);
			
			/* Promtly skip the event entirely in case a SCI cut isn't met. */
			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(sci21.hits.size() >= 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31_cut->Fill(sci31.E);

			const size_t N = foot->track.size();

			/* In this case, don't cut on any charges etc,.. just take the whole event and try to do
			 * general vertexing, angles, rho, etc. */
			if(!vertex_dist_cut_given) {
				h1_track_mult->Fill( foot->HasData()? N: -1 );

				std::vector<Line3D> tracks {};
				
				/* Sometime no selection is provided, in that case just fetch everything. */
				for(size_t i=0; i<N; ++i) {
					const RNFOOTTrack& t = foot->track[i];

					h2_q_vs_mult->Fill(N, t.Q);
					h2_score_vs_mult->Fill(N, t.score);

					tracks.push_back(*t);
				}
				if(N < 2) continue;

				const Point3D vertex = FindVertex(tracks);
				for(size_t i=0; i<N; ++i) {
					h2_track_distance->Fill( i, tracks[i].DistanceTo(vertex) );
				}
				double sum2 = 0;
				double theta;
				const Line3D heavy_track = tracks.front();
				for(size_t i=1; i < N; ++i) {
					theta = 1000.0 * tracks[i].AngleRelativeTo( heavy_track );
					sum2 += theta*theta;
					h1_track_angle->Fill( theta );
				}
				const double invariant_theta = sqrt(sum2);
				h1_angle_ex->Fill( invariant_theta );
				h2_vertex_z->Fill( vertex.z, invariant_theta);
			}

			else {
				VertexingResult<RNFOOTTrack> vtr =
					FindVertexingTracksMut(foot->track, distance_cut);

				if(!vtr.valid()) continue;

				const std::vector<RNFOOTTrack>& ftracks = vtr.tracks;
				const Point3D& vertex = vtr.vertex;
				const u32 n_tracks_selected = vtr.tracks.size();
				
				if(sum_n_tracks_required > 0 and (n_tracks_selected != sum_n_tracks_required))
					continue;
				
				u32 mask = (1U << n_tracks_selected) - 1; /* sequence 0b 1111.. */
				u32 unique_track_index = 0;
				/* Next, we require that this vertexing sequence of tracks also
				 * satisfies the charge cut. Selection windows are already
				 * presorted in descending charge. */
				for(const auto& [qcut, n_tracks_required] : selected) {
					A2 const& charge_cut = std::get<A2>( qcut );
					
					for(u32 i=0; i<n_tracks_required; ++i) {
						if( mnd::IsInside(
							ftracks.at(unique_track_index).Q,
							charge_cut))
						{ // Reset the unique bit.
							mask &= ~(1U << unique_track_index);
						}
						++unique_track_index;
					}
				}

				#ifdef MND_PHYSICS_PARANOIA
				assert( unique_track_index == sum_n_tracks_required
					&& "Huh seriously am I so bad?" );
				#endif
				/* `mask` now must be 0, otherwise we didn't catch all unique tracks for this
				 * specific charge interval. */
				if(sum_n_tracks_required > 0 && mask != 0)
					continue;
			
				/* Due to lexicographical combinations, the cute fact is that
				 * vtr.tracks is also already sorted accordingly :) */
				#ifdef MND_PHYSICS_PARANOIA
				assert( std::is_sorted(ftracks.begin(), ftracks.end())
					&& "Huh really?" );
				#endif
				
				/* At this point, selection is completed.
				 * Manifest back the Line3D objects to do angle calculations. */
				std::vector<Line3D> tracks;
				tracks.reserve(n_tracks_selected);
				for(const auto& ft : ftracks) {
					h2_q_vs_mult->Fill(n_tracks_selected, ft.Q);
					h2_score_vs_mult->Fill(n_tracks_selected, ft.score);
					tracks.push_back(*ft);
				}

				double sum2 = 0;
				double theta;
				switch(angle_type) {
					/* In this case, look for angles between
					* 'light' particles and the heaviest one. */
					case(AngleType::p): {
						const mnd::geom::Line3D heavy_track = tracks.front();
						for(u32 i=1; i < n_tracks_selected; ++i) {
							theta = 1000.0 * tracks[i].AngleRelativeTo( heavy_track );
							sum2 += theta*theta;
							h1_track_angle->Fill( theta );
						}
						break;
					}
					/* In this case, look for angles between
					 * 'heavy' particles and the lightest one. */
					case(AngleType::he): {
						const mnd::geom::Line3D light_track = tracks.back();
						for(u32 i=0; i < n_tracks_selected-1; ++i) {
							theta = 1000.0 * tracks[i].AngleRelativeTo( light_track );
							sum2 += theta*theta;
						}
						break;
					}
					/* In this case, angular `rho` is simply the RMS of
					 * all the possible angles.  */
					case(AngleType::all): {
						for(u32 i=0; i<n_tracks_selected; ++i) {
							for(u32 j=i+1; j<n_tracks_selected; ++j) {
								theta = 1000.0 * tracks[i].AngleRelativeTo( tracks[j] );
								sum2 += theta*theta;
								h1_track_angle->Fill( theta );
							}
						}
						break;
					}
				}

				const double invariant_theta = sqrt(sum2);
				h1_angle_ex->Fill( invariant_theta );
				h2_vertex_z->Fill( vertex.z, invariant_theta);

				for(u32 i=0; i<sum_n_tracks_required; ++i) {
					h2_track_distance->Fill( i, tracks[i].DistanceTo(vertex) );
				}

				h1_track_mult->Fill(n_tracks_selected);

			} // if(vertex_dist_cut_given)

			h1_sci21_cut2->Fill(sci21.E);
			h1_sci22_cut2->Fill(sci22.E);
			h1_sci31_cut2->Fill(sci31.E);

		} // for(size_t entryId{0}; entryId < nentries; ++entryId )

		bar.mark_as_completed();

	} // for(size_t i{0}; i < fileName.size(); ++i)

	show_console_cursor(true);

	TCanvas* cm = new TCanvas("Multp", "Recognized tracks", 2150, 1400);
	cm->Divide(2,2);
	cm->cd(1); h2_q_vs_mult->Draw("COLZ"); gPad->SetLogz();
	cm->cd(3); h2_score_vs_mult->Draw("COLZ"); gPad->SetLogz();
	cm->cd(4); h1_track_mult->Draw();
	cm->cd(2); new PLatex(0.08,
		"Coefficients: ",
		Form("Cr = %.1f mm^-2", Cr),
		Form("Cq = %.1f Q^-2", Cq),
		Form("Ct = %.1f mm^-2", Ct),
		Form("max cost for candidate: %.1f", max_cost),
		Form("max cost for whole track: %.1f", max_cost_f)
	);

	TCanvas* ct = new TCanvas("TrackDistance", "Recognized tracks distances between each other", 2150, 1400);
	ct->Divide(2,2);
	ct->cd(1); h2_track_distance->Draw();
	ct->cd(2); h1_track_angle->Draw();
	ct->cd(3); h1_angle_ex->Draw();
	ct->cd(4); h2_vertex_z->Draw("COLZ");

	TCanvas* cs = new TCanvas("SCIs", "SCI21,22,31", 2150, 1400);
	cs->Divide(3,3);
	cs->cd(1); h1_sci21->Draw();
	cs->cd(2); h1_sci22->Draw();
	cs->cd(3); h1_sci31->Draw();
	cs->cd(4); h1_sci21_cut->Draw();
	cs->cd(5); h1_sci22_cut->Draw();
	cs->cd(6); h1_sci31_cut->Draw();
	cs->cd(7); h1_sci21_cut2->Draw();
	cs->cd(8); h1_sci22_cut2->Draw();
	cs->cd(9); h1_sci31_cut2->Draw();

	WARN("Info: "); std::cerr << info << std::endl;

	canvas::save_all<canvas::Exe>( save, mnd::to_views(info) );

	WARN("End-of-main\n");
	rootApp.Run(); return 0;
}

std::istream& operator>>(std::istream& is, SingleSelect& out) {
	std::string s;
	is >> s;
	auto parts = mnd::split_view(s, ','); // vector<string_view>

	switch(parts.size()) { // {elem, n}
		case(2) : {
			auto z = magic_enum::enum_cast<AtomicNumber>(parts[0]);
			if(!z) {
				is.setstate(std::ios::failbit);
				return is;
			}

			u32 n{};
			if(!mnd::parse(parts[1], n)) {
				is.setstate(std::ios::failbit);
				return is;
			}

			out = SingleSelect {
				.fragment = *z,
				.index = n
			};
			break;
		}
		case(3) : {
			double a, b;
			u32 n;

			if(!mnd::parse(parts[0], a) ||
			   !mnd::parse(parts[1], b) ||
			   !mnd::parse(parts[2], n)) {
				is.setstate(std::ios::failbit);
				return is;
			}

			out = SingleSelect {
				.fragment = A2{a, b},
				.index = n
			};
			break;
		}

		default:
			is.setstate(std::ios::failbit);
	}
	return is;
}

std::ostream& operator<<(std::ostream& os, const SingleSelect& x) {
	switch(x.fragment.index()) {
		case(0):
			os << magic_enum::enum_name( std::get<0>(x.fragment) ); break;
		case(1):
			os << std::get<1>(x.fragment); break;
		default:
			__builtin_unreachable();
	}
	return os << ':' << x.index;
}
