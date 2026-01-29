#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

constexpr int xs[] = {1,3,5,7};
constexpr int ys[] = {0,2,4,6};

void foot_spread(std::string fileName = "", std::string o="x") {
	if(o != "x" and o != "y")
		throw std::runtime_error("Second argument must be either \"x\" or \"y\".");
	
	const int* Is = ((o == "x") ? &xs[0] : &ys[0]);

	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	std::array <
		std::shared_ptr<RNFOOTCal>, 4
	> foot {};

	for(int _i=0; _i<4; ++_i) {
		int i = Is[_i];
		foot[_i] = model->MakeField<RNFOOTCal>(Form("FOOT%d", i));
	}
	auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);

	TH1I *h2_foot[4];
	for(int i=0; i<4; ++i) {
		h2_foot[i] = new TH1I(Form("FOOT%d", Is[i]), Form("FOOT%d measuring %s", Is[i], o.c_str()), 320,0,640);
		h2_foot[i]->GetXaxis()->SetTitle(Form("%s in strip #", o.c_str()));
	}

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		for(int i=0; i<4; ++i) {
			std::vector<double> cl_pos = foot[i]->X();
			
			for(const auto p : cl_pos) {
				h2_foot[i]->Fill(p);
			}
		}
	}

	double fitMin = ((o == "x") ? 320 - 100 :  65);
	double fitMax = ((o == "x") ? 320 + 100 : 575);
	/* Try to fit a Gauss around them. */
	TCanvas* c = new TCanvas("c", "c", 1800, 1400);
	c->Divide(2,2);
	for(int i=0; i<4; ++i) {
		c->cd(i+1);
		gPad->SetGrid();

		TH1I* foot = h2_foot[i];
#define GAUS_PLUS_LIN
#ifdef GAUS_PLUS_LIN
		TF1 *f = new TF1(Form("f%d", i), "pol0(0) + gaus(1)", fitMin, fitMax);

		f->SetParameter(0, foot->GetMinimum());   // offset
		f->SetParameter(1, foot->GetMaximum());   // amplitude
		f->SetParameter(2, foot->GetMean());      // mean
		f->SetParameter(3, foot->GetRMS());       // sigma
		f->SetParLimits(3, 1e-6, 300);            // keep sigma positive
#else
		TF1 *f = new TF1(Form("f%d", i), "gaus(0) + gaus(1)", fitMin, fitMax);
		f->SetParameter(0, 0.75*foot->GetMaximum());  // amplitude
		f->SetParameter(3, 0.25*foot->GetMaximum());  // amplitude
		f->SetParameter(1, foot->GetMean());          // mean
		f->SetParameter(4, foot->GetMean());          // mean
		f->SetParameter(2, 0.9*foot->GetRMS());       // sigma
		f->SetParameter(5, 1.8*foot->GetRMS());       // sigma
#endif
		foot->Fit(f, "QR");
#ifdef GAUS_PLUS_LIN
		printf("FOOT%d : sigma = %.2f ± %.2f strips.\n", Is[i], f->GetParameter(3), f->GetParError(3));
#else
		printf("FOOT%d : sigma = %.2f ± %.2f strips.\n", Is[i], f->GetParameter(2), f->GetParError(2));
#endif

		foot->Draw("COLZ");
	}
}
