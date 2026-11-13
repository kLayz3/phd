/* This include is just to have LSP support when the define is present at compile-time
 * but not during editing. */
//#define MND_DEBUG_ENABLED

#include "util/FTrack.h"
#include "util/CLI.h"
#include "util/Geometry.h"
#include "util/PrettyHisto.h"
#include "util/MacroHelpers.h"
#include "util/TextUtils.h"

#include "TApplication.h"
#include "TParameter.h"
#include "TFOOTHitCont.h"
#include "TFRSHitCont.h"

using namespace ROOT;
using namespace ROOT::Experimental;
using namespace indicators;
using namespace mnd::col::literals;

int main(int argc, char* argv[]) {
	TClass* _cl = TClass::GetClass(typeid(RNFOOTTrack));
	if(!_cl or !_cl->GetDataMember("_u"))
		ERROR("MND_FOOTTRACK_DEBUG not compiled in, when the ROOT file got generated. Can't proceed\n");

	ROOT::EnableThreadSafety();
	CLI::App app{"Gate on straight primary beam tracks, and try to estimate the corresponding Cq, Cp, Ct coefficients, \
	              Based on already preliminary results from the Kalman filter. Basically, here we can gate explicitly \
                  on specific \'Z\' tracks (can be wrong z tho)"};
	
	std::vector<std::string> fileNames {};
	std::map<u32, A2> cut_q = {};
	A3 binning_q = {120, 0.0, 7.5};
	A3 binning_dx = {100,-5,5};
	A3 binning_dy = {100,-5,5};
	A3 binning_dq = {120, -2.1, 2.1};
	double fit_side_ratio = 2.0;
	A2 sci21_cut = {NAN,NAN};
	A2 sci22_cut = {NAN,NAN};
	A2 sci31_cut = {NAN,NAN};
	bool draw_sci = false;
	A3 sci_binning = {500, 300, 4000};

	std::vector<canvas::Extension> save = {};
	u32 nthreads = 1;
	
	add_logged_option(app, "-f,--file", fileNames, "Pass one or more file names, delimited by ','")
		->delimiter(',')
		->check(CLI::ReadPermissions);
	add_logged_option(app, "-q,--cut-q", cut_q, Form("Charge (Z) cut on each FOOT index, index=0,1,..%u.", N_FOOT_DETECTORS));
	add_logged_option(app, "-s,--side-ratio", fit_side_ratio, "When fitting the residues, select how much of the σ-width of "
		" the histogram should go into the fit. This is to enable precise targetted peak fitting, regardless of outlying bins.")
		->check(CLI::PositiveNumber);
	add_logged_option(app, "--binning-q", binning_q, "Binning of the charge measurement, Q of each FOOT, in [Z].")
		->delimiter(',');
	add_logged_option(app, "--binning-dx", binning_dx, "Binning of the difference δu (foot - track), in [mm].")
		->delimiter(',');
	add_logged_option(app, "--binning-dy", binning_dy, "Binning of the difference δv (foot - track), in [mm].")
		->delimiter(',');
	add_logged_option(app, "--binning-dq", binning_dq, "Binning of the difference δQ (foot - track), in [Z].")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci21", sci21_cut, "SCI21 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci22", sci22_cut, "SCI22 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option<DisplayDefault::No>(app, "--sci31", sci31_cut, "SCI31 QDC cut (also implying multiplicity 1). Default no cut.")
		->delimiter(',');
	add_logged_option(app, "--sci-binning",sci_binning, "Scintillator QDC binning.")
		->delimiter(',');
	add_logged_flag(app, "--draw-sci", draw_sci, "Draw sci histograms.");
	add_logged_option(app, "-o,--save", save, "Save the resulting histogram as an extension.")
		->delimiter(',');
	add_logged_option(app, "-n,--nthreads", nthreads,
		"Run the process multithreaded. One file per thread."
	)->check(CLI::PositiveNumber);
	
	bool test = false;
	add_logged_flag(app, "--test", test, "Test the CLI. Once parsed, just exit the program.");

	CLI11_PARSE(app, argc, argv);
	
	if(test) return 0;

#ifdef MND_FOOTTRACK_DEBUG
	if(fileNames.empty()) {
		WARN("To continue, must supply atleast one valid file name!\n"); return 0;
	}
	constexpr u32 N_PAIRS = RNFOOTHit::N_PAIRS;
	static_assert(N_PAIRS > 2, "To properly run this program must have N_PAIRS > 2.");
	
	double Cr, Cq, Ct, max_cost, max_cost_f;
	std::array<TH1I*, N_PAIRS> h1_diff_q;
	std::array<TH1I*, N_PAIRS> h1_diff_r;
	std::array<TH1I*, N_PAIRS> h1_diff_t;
	std::array<TH2I*, N_PAIRS> h2_diff_q;
	std::array<TH2I*, N_PAIRS> h2_diff_r;
	std::array<TH2I*, N_PAIRS> h2_diff_t;
	std::shared_ptr<TrackingGeometry<N_PAIRS>> geometry;
	{
		A3* c;
		TParameter<double>* m;
		auto f = std::make_unique<TFile>(fileNames.front().c_str(), "READ");
		get_obj(f, c, "FOOT_cost_coeff");
		get_obj(f, m, "FOOT_max_cost");
		Cr = c->at(0); Cq = c->at(1); Ct = c->at(2);
		max_cost = m->GetVal();
		get_obj(f, m, "FOOT_max_cost_f");
		max_cost_f = m->GetVal();
		for(u32 n=0; n<N_PAIRS; ++n) {
            get_obj(f, h1_diff_r.at(n), Form("FOOT_h1_diff_r_%u", n));
			get_obj(f, h1_diff_q.at(n), Form("FOOT_h1_diff_q_%u", n));
            get_obj(f, h1_diff_t.at(n), Form("FOOT_h1_diff_t_%u", n));

            get_obj(f, h2_diff_r.at(n), Form("FOOT_h2_diff_r_%u", n));
			get_obj(f, h2_diff_q.at(n), Form("FOOT_h2_diff_q_%u", n));
            get_obj(f, h2_diff_t.at(n), Form("FOOT_h2_diff_t_%u", n));
		}
		STrackingGeometry* stg;
		get_obj(f, stg, "FOOT_geometry");
		geometry = std::move( tg_from_stg(*stg) );
	}
	// Add all the other objects into h1_* and h2_* histograms. The cost values should be the same, regardless.
	for(size_t i=1; i<fileNames.size(); ++i) {
		
		auto f = std::make_unique<TFile>(fileNames[i].c_str(), "READ");
		for(u32 n=0; n<N_PAIRS; ++n) {
            auto h1_diff_r_tmp = get_obj<TH1I>(f, Form("FOOT_h1_diff_r_%u", n));
			auto h1_diff_q_tmp = get_obj<TH1I>(f, Form("FOOT_h1_diff_q_%u", n));
            auto h1_diff_t_tmp = get_obj<TH1I>(f, Form("FOOT_h1_diff_t_%u", n));

            auto h2_diff_r_tmp = get_obj<TH2I>(f, Form("FOOT_h2_diff_r_%u", n));
			auto h2_diff_q_tmp = get_obj<TH2I>(f, Form("FOOT_h2_diff_q_%u", n));
            auto h2_diff_t_tmp = get_obj<TH2I>(f, Form("FOOT_h2_diff_t_%u", n));

			h1_diff_r.at(n)->Add( h1_diff_r_tmp.get() );
			h2_diff_r.at(n)->Add( h2_diff_r_tmp.get() );
			h1_diff_q.at(n)->Add( h1_diff_q_tmp.get() );
			h2_diff_q.at(n)->Add( h2_diff_q_tmp.get() );
			h1_diff_t.at(n)->Add( h1_diff_t_tmp.get() );
			h2_diff_t.at(n)->Add( h2_diff_t_tmp.get() );
		}
	}
	
	std::vector<TH1P> resx, resy, resq;
	resx.reserve(N_PAIRS); resy.reserve(N_PAIRS); resq.reserve(N_PAIRS);
	
	std::vector<TH1P> h1_q, h1_dq;
	std::vector<TH2P> h2_qxy;
	h1_q.reserve(N_PAIRS); h1_dq.reserve(N_PAIRS);
	h2_qxy.reserve(N_PAIRS);
	std::vector<TH1P> h1_footx, h1_footy;
	h1_footx.reserve(N_PAIRS); h1_footy.reserve(N_PAIRS);
	for(u32 i=0; i<N_PAIRS; ++i) {
		resx.emplace_back(Form("((h1_fitrx_%d))Fit residue [mm]@FOOT%d u", i, i), kMagenta-9, binning_dx);
		resy.emplace_back(Form("((h1_fitry_%d))Fit residue [mm]@FOOT%d v", i, i), kYellow-9, binning_dy);
		resq.emplace_back(Form("((h1_fitrq_%d))Fit residue [charge unit]@FOOT%d Q", i, i), kCyan-9, binning_dq);

		h1_q.emplace_back(Form("((h1_q_%d))Q value@FOOT%d", i, i), 0xE3FF50_c, binning_q);
		h1_dq.emplace_back(Form("((h1_dq_%d))Difference Qx-Qy@FOOT%d", i, i), kBlue-9, binning_dq);
		h2_qxy.emplace_back(Form("((h2_q_%d))Q value y-detector [e]:Q value x-detector [e]@FOOT%d", i, i), binning_q, binning_q);

		h1_footx.emplace_back(Form("((h1_footx_%d))FOOT%d u [mm]", i, i), kMagenta-9, 320, -50, 50);
		h1_footy.emplace_back(Form("((h1_footy_%d))FOOT%d v [mm]", i, i), kMagenta-9, 320, -50, 50);
	}

	auto h2_kq = TH2P{"Kq cost [Q^2]:FOOT ID", N_PAIRS, -0.5, N_PAIRS-0.5, 500, 0, 1};
	auto h2_kr = TH2P{"Kr cost [mm^2]:FOOT ID", N_PAIRS, -0.5, N_PAIRS-0.5, 400, 0, 10};
	auto h2_kq_sqrt = TH2P{"sqrt_Kq cost [Q]:FOOT ID", N_PAIRS, -0.5, N_PAIRS-0.5, 500, 0, 1};
	auto h2_kr_sqrt = TH2P{"sqrt_Kr cost [mm]:FOOT ID", N_PAIRS, -0.5, N_PAIRS-0.5, 400, 0, 3.1628};
	auto h1_kt = TH1P{"Kt cost [mm^2]", kYellow-3, 400, 0, 10};
	auto h1_kt_sqrt = TH1P{"sqrt_Kt cost [mm]", kYellow-3, 400, 0, 3.162};

	auto h2_ktxy = TH2P{"FOOT - Upstream Y [mm]:FOOT - Upstream x [mm]@target", 400,-10,10, 400,-10,10};
	auto h1_diff_upstr_down = TH1P{"Distance Upstream to Downstream track [mm]", kGreen-1, 400, 0, 100};
	auto h1_has_upstream = TH1P{"Upstream track present [0=false, 1=true]", kYellow-3, 4, -0.75, 1.25};
	auto h1_score = TH1P{"Track score [a.u.]", kGreen-1, 500, 0, 50};

	auto h1_sci21 = TH1P{"SCI21 QDC mean [QDC units]", 0xCB00CB_c, sci_binning};
	auto h1_sci22 = TH1P{"SCI22 QDC mean [QDC units]", 0x0070DD_c, sci_binning};
	auto h1_sci31 = TH1P{"SCI31 QDC mean [QDC units]", 0xCC70DD_c, sci_binning};
	auto h1_sci21_cut  = TH1P{"((h1_cut)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, sci_binning};
	auto h1_sci22_cut  = TH1P{"((h1_cut)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, sci_binning};
	auto h1_sci31_cut  = TH1P{"((h1_cut)) SCI31 QDC mean [QDC units]@With cut", 0xCC80FD_c, sci_binning};
	auto h1_sci21_cut2 = TH1P{"((h1_cut2)) SCI21 QDC mean [QDC units]@With cut", 0x890389_c, sci_binning};
	auto h1_sci22_cut2 = TH1P{"((h1_cut2)) SCI22 QDC mean [QDC units]@With cut", 0x6180FD_c, sci_binning};
	auto h1_sci31_cut2 = TH1P{"((h1_cut2)) SCI31 QDC mean [QDC units]@With cut", 0xCC80FD_c, sci_binning};

	for(const auto& [foot_i, _] : cut_q) {
		if(foot_i >= N_PAIRS)
			ERROR("Supplied index: %u >= %u as Q-cut pair index.\n", foot_i, N_PAIRS);
	}
	
	TApplication rootApp("app", 0, 0);
	DynamicProgress<ProgressBar> progress;
	mnd::parallel_process(fileNames, nthreads, [=, &progress](size_t i, auto fname) mutable {
		auto model = RNTupleModel::Create();
		auto foot = model->MakeField<RNFOOTHit>("FOOT");
		auto frs = model->MakeField<RNFRSHit>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h104", fname);
		const size_t nentries = ntuple->GetNEntries();

		/* Fitting containers, for the track made by the other FOOT pairs. */
		std::array<double, N_PAIRS-1> qf;

		auto _bar_ptr = std::make_unique<ProgressBar>(
			option::BarWidth{50},
				option::Start{"["},
				option::Fill{"="},
				option::Lead{">"},
				option::Remainder{" "},
				option::End{"]"},
                option::PostfixText{mnd::msg(
					"%zu/%zu: %'zu (%s)", i+1, fileNames.size(), nentries, fname.c_str()
				)},
				option::ForegroundColor{ indicators::next_col() },
				option::ShowPercentage{true},
				option::ShowElapsedTime{true},
				option::ShowRemainingTime{true},
				option::FontStyles{std::vector{FontStyle::bold}}
		);
		auto idx = progress.push_back(std::move(_bar_ptr));
		auto& bar = progress[idx];

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
            if(mnd::PrintProgress(bar, entryId, nentries, 1000))
				progress.print_progress();

			const auto& sci21 = frs->cal.sci[RNFRSCal::SCI21_I];
			const auto& sci22 = frs->cal.sci[RNFRSCal::SCI22_I];
			const auto& sci31 = frs->cal.sci[RNFRSCal::SCI31_I];
			
			if(sci21.hits.size() >= 1) h1_sci21->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31->Fill(sci31.E);
			
			/* Skip the event entirely in case a SCI cut isn't met. */
			if(mnd::IsValid(sci21_cut) and (sci21.hits.size() != 1 or !mnd::IsInside(sci21.E, sci21_cut))) continue;
			if(mnd::IsValid(sci22_cut) and (sci22.hits.size() != 1 or !mnd::IsInside(sci22.E, sci22_cut))) continue;
			if(mnd::IsValid(sci31_cut) and (sci31.hits.size() != 1 or !mnd::IsInside(sci31.E, sci31_cut))) continue;

			if(sci21.hits.size() >= 1) h1_sci21_cut->Fill(sci21.E);
			if(sci22.hits.size() >= 1) h1_sci22_cut->Fill(sci22.E);
			if(sci31.hits.size() >= 1) h1_sci31_cut->Fill(sci31.E);

			const double x0_upst = frs->xT;
			const double y0_upst = frs->yT;
			
			const mnd::geom::Line3D lu = RNTrackToLine3D( frs->s2_bt );
			h1_has_upstream->Fill( static_cast<int>(lu.IsValid()) );

			if(foot->track.empty()) continue;
			
			bool has_valid_track = false;

			/* Find the corresponding FOOT track(s). */
			for(const auto& track: foot->track) {
				if(track.n != N_PAIRS) continue;

				/* Check if charge veto passed. Must pass for all the detectors. */
				if(!std::all_of(cut_q.begin(), cut_q.end(),
					[&track](const auto& e) -> bool {
						const auto& [foot_i, cut] = e;
						const double qx = track._qx[ foot_i ];
						const double qy = track._qy[ foot_i ];
						const double q  = (qx + qy)/2;
						return mnd::IsInside(q, cut);
					})
				) {
					continue;
				}
				has_valid_track = true;

				// Loop over individual layers, check the residues..
				for(u32 i=0; i<N_PAIRS; ++i) {
					const double qx = track._qx[i];
					const double qy = track._qy[i];
					const double qmean = (qx + qy) / 2;
					const double qdiff = qx - qy;
					const double qvar = 0.5*qdiff*qdiff;

					// Mask for which dets fired is all except one at bit i.
					const u32 mask = ((1 << N_PAIRS) - 1) & ~(1 << i);
					
					const auto& parameters =
						geometry->fit.get<N_PAIRS - 1>().get_from_mask(mask);
					using Scalar = FLayerInfo::EntryXYUnderlyingType; // double|long double
					Eigen::Vector<Scalar, 2 * (N_PAIRS - 1)> rhs;
					u32 w = 0;

					/* Populate the arrays for N-1 fit. */
					for(u32 j=0; j<N_PAIRS; ++j) {
						if(j == i) continue;

						rhs(2*w) = track._u[j];
						rhs(2*w + 1) = track._v[j];
						qf[w] = (track._qx[j] + track._qy[j]) / 2;
						++w;
					} if(w != N_PAIRS-1) ERROR("w must be %u, hello?", N_PAIRS-1);
					
					/* For the linear design problem, e.g. measurements to (x0,ax,y0,ay) track
					* need to add the offsets to the vector side. */
					rhs += parameters.second;

					const Eigen::Vector<Scalar, 4> normal = parameters.first.solve(rhs);
					const Eigen::Vector<Scalar, 2> predicted =
						geometry->layer[i].A * normal - geometry->layer[i].dxy;
					auto [q_extr, track_qvar] = mnd::mean_var(qf);

					double dx = track._u[i] - predicted.x();
					double dy = track._v[i] - predicted.y();
					double dq = qmean - q_extr;

					h1_footx[i].Fill( track._u[i] );
					h1_footy[i].Fill( track._v[i] );
					resx[i].Fill(dx);
					resy[i].Fill(dy);
					resq[i].Fill(dq);

					h1_q[i].Fill(qmean);
					h1_dq[i].Fill(qdiff);
					h2_qxy[i].Fill(qx, qy);
					
					/* Find the 3-measurement 'test' track charge params.. */
					h2_kq.Fill( i, dq*dq + qvar );
					h2_kr.Fill( i, dx*dx + dy*dy );
					h2_kq_sqrt.Fill( i, sqrt(dq*dq + qvar) );
					h2_kr_sqrt.Fill( i, sqrt(dx*dx + dy*dy) );
				} // for(u32 i=0; i<N_PAIRS; ++i)

				const mnd::geom::Line3D ld = RNTrackToLine3D(track);
				const mnd::geom::Point2D r0_down = ld.Eval( TFOOTHitCont::TARGET_Z );

				const double dx_dstr_ustr = r0_down.x - x0_upst;
				const double dy_dstr_ustr = r0_down.y - y0_upst;
				
				if(lu.IsValid() && std::isfinite(dx_dstr_ustr) && std::isfinite(dy_dstr_ustr)) {
					h2_ktxy->Fill(dx_dstr_ustr, dy_dstr_ustr);
					h1_kt->Fill(dx_dstr_ustr*dx_dstr_ustr + dy_dstr_ustr*dy_dstr_ustr);
					h1_kt_sqrt->Fill(sqrt(dx_dstr_ustr*dx_dstr_ustr + dy_dstr_ustr*dy_dstr_ustr));
					h1_diff_upstr_down->Fill(ld.DistanceTo(lu));
				}
				h1_score->Fill(track.score);

			}
			if(has_valid_track) {
				h1_sci21_cut2->Fill(sci21.E);
				h1_sci22_cut2->Fill(sci22.E);
				h1_sci31_cut2->Fill(sci31.E);
			}
		}

		if(!bar.is_completed())
			bar.mark_as_completed();
    });

	std::vector<std::string> resolutions {};
	
	/* Around each residue, also fit a teeny-weeny gauss-chan 🥺 👉👈
	 * Its' width is what we label as detector resolution. */
	TCanvas* c1d = new TCanvas("FitResidue1D", "Fit residues 1D", 2150, 1650);
	c1d->Divide(N_PAIRS, 3);
	for(size_t i=0; i<N_PAIRS; ++i) {
		c1d->cd(i+1);
		auto [fit_result_x, _x_] = resx[i].DrawAndFit(fit_side_ratio, kMagenta + 2, 5, 4 /* niter */);
		WARN("FOOT%zu X: " KBH_MAG "%.2f, %.2f" KNRM " [um]\n", i, 1e3*fit_result_x[1], 1e3*fit_result_x[2]);
		c1d->cd(i+1 + N_PAIRS);
		auto [fit_result_y, _y_] = resy[i].DrawAndFit(fit_side_ratio, kMagenta + 2, 5, 4 /* niter */);
		WARN("FOOT%zu Y: " KBH_YEL "%.2f, %.2f" KNRM " [um]\n", i, 1e3*fit_result_y[1], 1e3*fit_result_y[2]);

		c1d->cd(i+1 + 2*N_PAIRS);
		auto [fit_result_q, _q_] = resq[i].DrawAndFit(fit_side_ratio, kMagenta + 2, 5, 4 /* niter */);
		WARN("FOOT%zu Q: " KBH_CYN "%.3f, %.3f" KRNM " [Q]\n", i, fit_result_q[1], fit_result_q[2]);

		resolutions.emplace_back (
			Form("FOOT%zu (X,Y): (%s%.2f%s, %s%.2f%s) um",
				i, KRED, 1e3*fit_result_x[2], KNRM, KBLU, 1e3*fit_result_y[2], KNRM)
		);
	}
	TCanvas* cpos = new TCanvas("foot_xyq", "FOOT absolute positions/charge", 2150, 1650);
	cpos->Divide(N_PAIRS, 2);
	for(size_t i=0; i<N_PAIRS; ++i) {
		cpos->cd(i+1+ 0*N_PAIRS); h1_footx[i].Draw();
		cpos->cd(i+1+ 1*N_PAIRS); h1_footy[i].Draw();
	}
	TCanvas *cq = new TCanvas("foot_q", "FOOT charge-related histograms", 2150, 1650);
	cq->Divide(N_PAIRS, 3);
	for(size_t i=0; i<N_PAIRS; ++i) {
		cq->cd(i+1+ 0*N_PAIRS); h1_q[i].Draw();
		cq->cd(i+1+ 1*N_PAIRS); h1_dq[i].Draw();
		cq->cd(i+1+ 2*N_PAIRS); h2_qxy[i].Draw("COLZ");
	}

	TCanvas* ckv = new TCanvas("CostValues2D", "Cost Values 2D", 2200, 1400);
	ckv->Divide(4, 2);
	ckv->cd(1); h2_kq.Draw("COLZ");
	ckv->cd(5); h2_kq_sqrt.Draw("COLZ");
	ckv->cd(2); h2_kr.Draw("COLZ");
	ckv->cd(6); h2_kr_sqrt.Draw("COLZ");
	ckv->cd(3); h1_kt.Draw();
	ckv->cd(7); h1_kt_sqrt.Draw();

	ckv->cd(4); h2_ktxy.Draw("COLZ");
	ckv->cd(8); h1_diff_upstr_down.Draw();

	TCanvas* cscore = new TCanvas("Score", "Score of the recognized tracks", 1200, 800);
	cscore->Divide(2,2);
	cscore->cd(1); h1_score.Draw();
	cscore->cd(2); new PLatex(0.08,
		"Coefficients: ",
		Form("Cr = %.1f", Cr),
		Form("Cq = %.1f", Cq),
		Form("Ct = %.1f", Ct),
		Form("max cost: %.1f", max_cost),
		Form("max cost_f: %.1f", max_cost_f)
	);
	cscore->cd(3); h1_has_upstream.Draw();
	cscore->cd(4); new PLatex(0.08,
		"Position resolution: ",
		mnd::as_span(resolutions)
	);

	TCanvas* cRQT = new TCanvas("R-Q-T-Score", "Individual R,Q,T-score components", 2150, 1400);
	cRQT->Divide(N_PAIRS, 3);
	TCanvas* cRQT2 = new TCanvas("R-Q-T-Score, 2D", "Individual R,Q,T-score components in 2D", 2150, 1400);
	cRQT2->Divide(N_PAIRS, 3);
#define DECORATE_HISTO(hshort, col) \
	for(u32 n=0; n<N_PAIRS; ++n) { \
		TH1I* h = h1_diff_##hshort.at(n); \
		h->SetFillColor(col); h->SetLineColor(kBlack); h->SetLineWidth(1);  \
	}
	DECORATE_HISTO(r, kRed-7);
	DECORATE_HISTO(q, kCyan-4);
	DECORATE_HISTO(t, kYellow-7);
	for(u32 n=0; n<N_PAIRS; ++n) {
		cRQT ->cd(0*N_PAIRS + 1+n); h1_diff_r.at(n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
		cRQT ->cd(1*N_PAIRS + 1+n); h1_diff_q.at(n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
		cRQT ->cd(2*N_PAIRS + 1+n); h1_diff_t.at(n)->Draw(); gPad->SetGrid(); gPad->SetLogy();
		cRQT2->cd(0*N_PAIRS + 1+n); h2_diff_r.at(n)->Draw("COLZ"); gPad->SetLogz();
		cRQT2->cd(1*N_PAIRS + 1+n); h2_diff_q.at(n)->Draw("COLZ"); gPad->SetLogz();
		cRQT2->cd(2*N_PAIRS + 1+n); h2_diff_t.at(n)->Draw("COLZ"); gPad->SetLogz();
	}

	if(draw_sci) {
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
	}

	double qvg = cut_q.size() > 0
	? (std::accumulate(cut_q.begin(), cut_q.end(), 0.0,
		[](double sum, const auto& value) -> double {
			return sum + (value.second[0] + value.second[1]) / 2;
		}) / cut_q.size() )
	: 0.0;

	canvas::save_all<canvas::Exe>(save, {
		mnd::rootfile_names_concatenated(fileNames),
		std::string{"Q:"} + mnd::ftos(qvg)
	});

	WARN("End-of-main");
	rootApp.Run();

#else
    WARN("\n" BOLD ".. MND_FOOTTRACK_DEBUG not compiled in. This program does nothing.\n" KNRM);
#endif
    return 0;
}
