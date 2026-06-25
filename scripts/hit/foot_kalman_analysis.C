#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"

#include "../../includes/util/Geometry.h"
#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/Tracking.h"
#include "../../includes/util/MacroHelpers.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

enum AtomicNumber { H, He, Li, Be, B, C };
using AllowOthers = DoSave; // also yes/no

using DistanceCut  = mnd::InputWrapper<A2>;
using FileSequence = mnd::InputWrapper<std::vector<std::string>>;
using FromFile     = mnd::InputWrapper<std::string>;
using Select = mnd::InputWrapper<
	std::vector<
		std::pair <
			std::variant<AtomicNumber, A2>, 
			uint32_t
		>
	>
>;

constexpr A2 elem_to_a2(AtomicNumber e) {
	switch(e) {
		case H:  return A2{0.6, 1.5};
		case He: return A2{1.5, 2.6};
		case Li: return A2{2.6, 3.5};
		case Be: return A2{3.5, 4.5};
		case B:  return A2{4.5, 5.4};
		case C:  return A2{5.4, 6.8};
	}
}	

void foot_kalman_analysis (
	std::variant<FileSequence, FromFile> fileNamesA = {},
	Select selected = {},
	AllowOthers allow_others = AllowOthers::yes,
	DistanceCut distance_cut = {0.0, 1000},
	DoSave do_save = DoSave::no
) {
	std::vector<std::string> fileNames;
	if( std::holds_alternative<FromFile>(fileNamesA) ) {
		fileNames = ParseFile( std::get<FromFile>(fileNamesA) );
	} else {
		fileNames = std::move( std::get<FileSequence>(fileNamesA) );
	}
	
	if(fileNames.size() == 0) return;
	WARN("Handling files: "); std::cerr << fileNames << std::endl;

	constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;
	double Cr, Cq, Ct, Cp, max_cost;
	{
		const auto& fileName = fileNames.front();
		std::array<double, 4>* c;
		TParameter<double>* m;
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		get_obj(f, c, "FOOT_cost_coeff");
		get_obj(f, m, "FOOT_max_cost");
		Cr = c->at(0); Cq = c->at(1); Ct = c->at(2); Cp = c->at(3);
		max_cost = m->GetVal();
	}
	/* Sanitize some cmd-line-args. */	
	A2 distance_cut_arr = static_cast<A2>(distance_cut);	
	u32 sum_n_tracks_required = std::accumulate(
		selected.begin(), selected.end(), 0, [](u32 sum, const auto& x) {
			return sum + x.second;
		}
	);
	if(sum_n_tracks_required == 0) 
		allow_others = AllowOthers::yes;

	/* Convert the selected enum into the A2. */
	for(auto& [key,value] : selected) {
		if(std::holds_alternative<AtomicNumber>(key)) {
			key = elem_to_a2( std::get<AtomicNumber>(key) );
		}
		if(value == 0) ERROR("Must require 1 or more tracks, cannot be 0. Stop that.");
	}
	// All cut entries now hold only A2 range.
	std::array<bool, 32> track_used_mask = {};

	TH1P* h1_track_mult = new TH1P("Track multiplicity", kRed-1, 10, -0.5, 9.5);
	TH2P* h2_q_vs_mult = new TH2P("Track charge [Q]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 40, 0,8);
	TH2P* h2_score_vs_mult = new TH2P("Track score [a.u.]:Track multp@Full FOOT system",
		10, -0.5, 9.5, 300, 0, 50);
	TH1P* h1_track_distance = new TH1P("Track distances [mm]", kBlue-1, 600, 0, 30);
	TH1P* h1_track_angle = new TH1P("Track angles [mrad]", kMagenta+1, 200, 0, 100);
	TH1P* h1_angle_ex = new TH1P("Track angles, excitation sqrt(t1^2 + t2^2 + t3^2) [mrad]", kCyan-9, 200, 0, 1000);
	TH2P* h2_vertex_z = new TH2P("RMS angle [mrad]:Vertex z [mm]@Traced by the FOOT", 200, -100, 100, 100,0,100);

	for(const auto& fname : fileNames) {
		WARN("Proceeding with file: \'%s\'\n", fname.c_str());
		auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTHit>("FOOT");
		auto frs = model->MakeField<RNFRSHit>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);

		ROOT::EnableImplicitMT();

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			bool is_valid = true;
			
			const size_t N = foot->track.size();
			if(N > track_used_mask.size()) continue; // abnormally large event?
			
			/* In case we don't allow rogue tracks, then the number of validated tracks 
			 * must be the total number. */
			if(allow_others == AllowOthers::no && N != sum_n_tracks_required) continue;
			
			/* Select the event only if corresponing charge values are found, in their
			 * respective 'quantity' */
			track_used_mask.fill(false);
			
			for(const auto& [cut, n_tracks_required] : selected) {
				const A2& charge_cut = std::get<A2>( cut ); 
				/* Go over the tracks, if it's sorted within a cut then 'mask' the correspodning entry for
				 * next iterations. */
				u32 n = n_tracks_required;
				for(size_t i=0; i<N; ++i) {
					if( track_used_mask[i] ) continue; // skip event. Entry already "sorted out" 
					const RNFOOTTrack& t = foot->track[i];
					if(mnd::IsInside(t.Q, charge_cut)) {
						track_used_mask[i] = true;
						if(--n == 0) break;
					};
				} 

				/* `n` now must be 0, otherwise we didn't catch all the required charges. */
				if(n != 0) is_valid = false;
			}
			if(!is_valid) continue; // skip event.
			
			h1_track_mult->Fill(N);

			/* From the selected tracks, fill out the `tracks` 
			 * sequence. Here, we eliminate tracks from the event that didn't fell in the 
			 * `selected` window(s). */
			std::vector<mnd::geom::Line3D> tracks {};
			
			/* Sometime no selection is provided, in that case just fetch everything. */
			for(size_t i=0; i<N; ++i) {
				/* If the current track wasn't marked by the `selected` filtering above, means it didn't
				 * fall into any of the selection windows. */
				if( selected.get().size() > 0 and !track_used_mask[i] ) continue; 
				const RNFOOTTrack& t = foot->track[i];
				
				h2_q_vs_mult->Fill(N, t.Q);
				h2_score_vs_mult->Fill(N, t.score);

				tracks.push_back( RNTrackToLine3D(t) );
			}
			const size_t N_valid = tracks.size();
		
			if(N_valid > 1) {
				/* Find the vertex. */
				const mnd::geom::Point3D vertex = mnd::geom::FindVertex(mnd::as_span(tracks));

				/* Cut on vertex distances... */
				std::vector<double> distances {};

				for(const mnd::geom::Line3D& track : tracks) {
					double distance = track.DistanceTo( vertex );
					if(mnd::IsValid(distance_cut) and !mnd::IsInside(distance, distance_cut)) {
						is_valid = false;
						break;
					}
					distances.push_back(distance);
				}
				// If a single track is not within with distance cut, skip event.
				if(!is_valid) continue;

				for(auto d: distances) h1_track_distance->Fill( d );

				double sum2 = 0;
				double sum = 0;
				for(size_t i=0; i<N_valid; ++i) {
					for(int j=i+1; j<N_valid; ++j) {
						double theta = 1000.0 * tracks[i].AngleRelativeTo( tracks[j] );
						sum2 += theta*theta;
						h1_track_angle->Fill( theta );
					}
				}
				h1_angle_ex->Fill( sqrt(sum2) );
				h2_vertex_z->Fill( vertex.z, sqrt(sum2) / N_valid );
			}
		}

		ROOT::DisableImplicitMT();
	}


	TCanvas* cm = new TCanvas("Multp", "Recognized tracks", 2150, 1400);
	cm->Divide(2,2);
	cm->cd(1); h2_q_vs_mult->Draw("COLZ"); gPad->SetLogz();
	cm->cd(3); h2_score_vs_mult->Draw("COLZ"); gPad->SetLogz();
	cm->cd(4); h1_track_mult->Draw();
	cm->cd(2); PLatex(0.08,
		"Coefficients: ",
		Form("Cr = %.1f mm^-2", Cr),
		Form("Cq = %.1f Q^-2", Cq),
		Form("Ct = %.1f mm^-2", Ct),
		Form("Cp = %.1f", Cp),
		Form("max cost [@4th]: %.1f", max_cost)
	);

	TCanvas* ct = new TCanvas("TrackDistance", "Recognized tracks distances between eachother", 2150, 1400);
	ct->Divide(2,2);
	ct->cd(1); h1_track_distance->Draw();
	ct->cd(2); h1_track_angle->Draw();
	ct->cd(3); h1_angle_ex->Draw();
	ct->cd(4); h2_vertex_z->Draw("COLZ");

	if(do_save == DoSave::yes) {
		std::stringstream pname;
		for(const auto& fname : fileNames) {
			std::filesystem::path p{fname};
			pname << p.replace_extension().c_str() << "--";
		}
		std::filesystem::path inf( pname.str() );
		save_all(canvas::Extension::png, { inf.stem().c_str() });
	}
}
