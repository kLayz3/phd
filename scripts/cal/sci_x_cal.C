#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

#include "../../includes/PolyFitter.hxx"
#include "../../includes/PrettyHisto.hxx"
constexpr const char* label[] = {
	"21", "22", "23", "24"
}; constexpr int NSCI = static_cast<int>( sizeof(label)/sizeof(*label) );

void sci_x_cal(std::string fileName = "", int i_sci = 0, double x_min = -DBL_MAX, double x_max = DBL_MAX) {
	if(i_sci != 0 and i_sci != 1)
		throw std::invalid_argument("Second arg (i_sci) must be either 0 or 1. Calibrating only SCI21 and SCI22 here.");

	ROOT::EnableImplicitMT();
	
	std::array<TPCParam, RNFRSCal::N_VALID_TPC> *tpc_params;
	std::array<SCIParam, RNFRSCal::N_VALID_SCI> *sci_params;
	{
		std::unique_ptr<TFile> f = std::make_unique<TFile>(fileName.c_str(), "READ");
		tpc_params = f->Get < 
			std::remove_reference_t<decltype(*tpc_params)>
		> ("FRS_tpc_parameters");
		if(!tpc_params)
			throw std::runtime_error(Form("TPC param is nullptr. Fix it (line: %d).", __LINE__));
		sci_params = f->Get < 
			std::remove_reference_t<decltype(*sci_params)>
		> ("FRS_sci_parameters");
		if(!sci_params)
			throw std::runtime_error(Form("SCI param is nullptr. Fix it (line: %d).", __LINE__));

	}
	const auto& sci_param = sci_params->at(i_sci); 
	const double z0 = sci_param.z0; 
	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSCal>("FRS");
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);


	/* Take TPC21,22,23,24 for the extrapolations. */
	constexpr int N = 4;
	const std::array<double, N> zTPC = {
		tpc_params->at(0).z0,
		tpc_params->at(1).z0,
		tpc_params->at(2).z0,
		tpc_params->at(3).z0,
	};
	/* Optionally skip some TPC's. */
//#define SKIP_TPC21
#define SKIP_TPC22
//#define SKIP_TPC23
#define SKIP_TPC24
	
	std::vector<double> xs, ys, zs;
	xs.reserve(N); ys.reserve(N); zs.reserve(N);
	std::array<double, 2> px, py;

	TH2P* hist  = new TH2P(Form("SCI%s X[mm]:TPC extr. [mm]", label[i_sci]), 200, -50, 50, 400, -100, 100); 
	TH2P* histd = new TH2P(Form("SCI%s X -  TPC extr.[mm]:TPC extr. [mm]", label[i_sci]), 200, -50, 50, 200, -10, 10); 
	TH1P* href  = new TH1P(Form("TPC extr. to SCI%s [mm]", label[i_sci]), ORGB{0xf0A330C9}, 200, -50, 50);
	TH2P* histy = new TH2P(Form("SCI%s X -  TPC extr.[mm]:TPC extr (Y). [mm]", label[i_sci]), 200, -50, 50, 200, -10, 10); 
	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		xs.clear(); ys.clear(); zs.clear();
		
		const auto& sci = frs->sci[i_sci];
		if(sci.hits.size() == 0) continue;

		for(int i = 0; i < N; ++i) {
#ifdef SKIP_TPC21
			if(i == 0) continue;
#endif
#ifdef SKIP_TPC22
			if(i == 1) continue;
#endif
#ifdef SKIP_TPC23
			if(i == 2) continue;
#endif
#ifdef SKIP_TPC24
			if(i == 3) continue;
#endif
			const auto& tpc = frs->tpc[i];
			if(tpc.hits[0].size() != 1 or tpc.hits[1].size() != 1) continue;

			double x = tpc.X0();
			double y = tpc.Y0();
			if(!std::isfinite(x) or !std::isfinite(y)) continue;
			
			// Limit x measurements of TPC's to a small value around 0. If supplied
			if(x < x_min or x > x_max) continue;

			zs.push_back( zTPC[i] );	
			xs.push_back(x);
			ys.push_back(y);
		}
		if(zs.size() < 2) continue;
		PolyFit<1>(zs, xs, px);
		PolyFit<1>(zs, ys, py);
		double x_extr = px[0] + px[1] * z0;
		double y_extr = py[0] + py[1] * z0;
		hist ->Fill(x_extr, sci.hits[0].x);
		histd->Fill(x_extr, sci.hits[0].x - x_extr);
		histy->Fill(y_extr, sci.hits[0].x - x_extr);
		href->Fill(x_extr);
	}

	TProfile* profile = (*hist)->ProfileX("_pf_x");	
	profile->Fit("pol1", "Q", "", -11, 11);
	TF1* fit = profile->GetFunction("pol1");

	if(fit) {
		double slope  = fit->GetParameter(1);
		double offset = fit->GetParameter(0);
		printf("Found: SCI%s: slope = %.5f [TDC_VAL/mm]; off = %.5f [TDC_VAL]. Inverse is: %.7f [mm/TDC_VAL], offset = %.7f [mm] \n", 
			label[i_sci], slope, offset, 1.0/slope, offset/slope);
		WARN("Recommended: " EBOLD(x_factor: %.5f) ", and " EBOLD(x_offset: %.5f\n),
				sci_param.x_factor / slope, (sci_param.x_offset - offset) / slope );
	} else {
		printf("Fit completely failed!\n");
	}
	printf("Currently in the file: (offset,slope) = (%.4f, %.4f)\n", sci_param.x_offset, sci_param.x_factor);
	TCanvas *c = new TCanvas("c", "c", 2400, 1400);
	c->Divide(2,2);
	c->cd(1); hist->Draw("COLZ");
	c->cd(2); histd->Draw("COLZ");
	c->cd(3); href->Draw();
	std::string s = std::string{"Extrapolation from: TPC"} +
#ifndef SKIP_TPC21
		"21 " +
#endif
#ifndef SKIP_TPC22
		"22 " +
#endif
#ifndef SKIP_TPC23
		"23 " +
#endif
#ifndef SKIP_TPC24
		"24 " +
#endif
		+ "";

	c->cd(4);
	PLatex(0.08, 
		std::move(s), 
#ifndef SKIP_TPC21
		Form("TPC21: %.1f", tpc_params->at(0).z0),
#endif
#ifndef SKIP_TPC22
		Form("TPC22: %.1f", tpc_params->at(1).z0),
#endif
#ifndef SKIP_TPC23
		Form("TPC23: %.1f", tpc_params->at(2).z0),
#endif
#ifndef SKIP_TPC24
		Form("TPC24: %.1f", tpc_params->at(3).z0),
#endif
		Form("SCI%s: %.1f", label[i_sci], sci_param.z0)
	);
	
	TCanvas *cy = new TCanvas("cy", "cy", 1200, 800);
	histy->Draw("COLZ");
}
