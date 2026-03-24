#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
using namespace ROOT;
using namespace ROOT::Experimental;

#include "../../includes/Plotter.hxx"

constexpr int N = RNFRSMap::N_VALID_SCI; 
const char* labels[N] = {
	"SCI21", "SCI22", "SCI31", "SCI41"	
};

void sci_cal_single(std::string, int );

void sci_cal(std::string fileName = "", int isci=-1) {
	if(isci >= N)
		throw std::invalid_argument(Form("Sci index: %d must be < %d", isci, N));
	if(isci < 0) {
		for(int i=0; i<N; ++i) {
			sci_cal_single(fileName, i);
		}
	}
	else
		sci_cal_single(fileName, isci);
}

void sci_cal_single(std::string fileName, int isci) {
	if(isci >= N or isci < 0)
		throw std::invalid_argument(Form("Second arg must be 0 <= x < %d index of sci.", N));
	
	ROOT::EnableImplicitMT();

	auto model = RNTupleModel::Create();
	auto frs = model->MakeField<RNFRSMap>("FRS"); // shared_ptr.
	auto ntuple = RNTupleReader::Open(std::move(model), "h102", fileName);

	TH1I* h1_ldiff = new TH1I(Form("h1_ldiff_%d", isci), Form("%s raw TDC left", labels[isci]), 1000,-10000, 200000);
	TH1I* h1_rdiff = new TH1I(Form("h1_rdiff_%d", isci), Form("%s raw TDC right", labels[isci]), 1000,-10000, 200000);
	TH1I* h1_lq = new TH1I(Form("h1_lq_%d", isci), Form("%s raw QDC left", labels[isci]), 4095,0,4095);
	TH1I* h1_rq = new TH1I(Form("h1_rq_%d", isci), Form("%s raw QDC right", labels[isci]), 4095,0, 4095);
	TH1I* h1_lr = new TH1I(Form("h1_lr_%d", isci), Form("%s raw TDC left-right", labels[isci]), 4000, -2000, 2000);
	TH2I* h2_qq_lr = new TH2I(Form("h1_qq_lr_%d", isci), Form("%s QDC (mean) vs. TDC diff left-right", labels[isci]), 4000, -2000, 2000, 4095,0,4095);
	
	h1_ldiff->GetXaxis()->SetTitle("Raw PMT left timing [25 ps]");
	h1_rdiff->GetXaxis()->SetTitle("Raw PMT left timing [25 ps]");
	h1_ldiff->GetYaxis()->SetTitle("Count");
	h1_rdiff->GetYaxis()->SetTitle("Count");
	h1_ldiff->SetFillStyle(1001); h1_ldiff->SetFillColor(kGreen + 1);
	h1_rdiff->SetFillStyle(1001); h1_rdiff->SetFillColor(kGreen - 1);

	h1_lq->GetXaxis()->SetTitle("Charge integral [ADC units]");
	h1_rq->GetXaxis()->SetTitle("Charge integral [ADC units]");
	h1_lq->GetYaxis()->SetTitle("Count");
	h1_rq->GetYaxis()->SetTitle("Count");
	h1_lq->SetFillStyle(1001); h1_lq->SetFillColor(kRed + 1);
	h1_rq->SetFillStyle(1001); h1_rq->SetFillColor(kRed - 1);

	h1_lr->GetXaxis()->SetTitle("TDC left - TDC right [25 ps]");
	h1_lr->GetYaxis()->SetTitle("Count");
	h1_lr->SetFillStyle(1001); h1_lr->SetFillColor(kCyan + 1);

	h2_qq_lr->GetXaxis()->SetTitle("TDC left - TDC right [25 ps]");
	h2_qq_lr->GetYaxis()->SetTitle("#surd(q_{1}#times q_{2}) [ADC units]");

	for(auto entryId : *ntuple) {
		ntuple->LoadEntry(entryId);
		
		const auto& sci = frs->sci[isci];
		const auto& tdc = sci.tdc;

		if(tdc.size() != 1) continue;

		int tdc_l = tdc[0].tdc_l;
		int tdc_r = tdc[0].tdc_r;

		h1_ldiff->Fill(tdc_l);
		h1_rdiff->Fill(tdc_r);
		h1_lq->Fill(sci.qdc[0]);
		h1_rq->Fill(sci.qdc[1]);
		h1_lr->Fill(tdc_l - tdc_r);
		h2_qq_lr->Fill(tdc_l - tdc_r , sqrt(sci.qdc[0]*sci.qdc[1]));
	}
	
	TCanvas* c = new TCanvas(Form("c%d", isci), Form("c%d", isci), 1800, 1250);
	c->Divide(3,2);
	//mnd::SetTitle(c, Form("map:%s", labels[isci]) + "-" + mnd::fname_short(fileName));
	mnd::SetTitle(c, mnd::sstrcat( Form("map:%s", labels[isci]), " - ", mnd::fname_short(fileName) ));

#define HILFE(h,i,...) \
	c->cd(i); \
	gPad->SetGrid(); \
	h->Draw(__VA_ARGS__);

	HILFE(h1_ldiff, 1, "HIST")
	HILFE(h1_rdiff, 2, "HIST")
	HILFE(h1_lq, 3, "HIST")
	HILFE(h1_rq, 4, "HIST")
	HILFE(h1_lr, 5, "HIST")
	HILFE(h2_qq_lr, 6, "COLZ")
}
