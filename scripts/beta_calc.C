#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../includes/PrettyHisto.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;

constexpr const char* label[] = {
	"21", "22", "23", "24"
}; constexpr int NSCI = static_cast<int>( sizeof(label)/sizeof(*label) );

constexpr uint32_t col[] = {
	0xAAD372, 0x3FC7EB, 0xF48CBA, 0xFFF468, 0x8788EE, 0xC69B6D
}; constexpr int NCOL = static_cast<int>( sizeof(label)/sizeof(*label) );

constexpr std::pair<int,int> comb[] = {
	{0,1}, /* ToF 21-22 */
	{1,2}, /* ToF 22-23 */
}; constexpr int NCOMB = static_cast<int>( sizeof(comb)/sizeof(*comb) );

constexpr int CANVAS_LENGTH = std::max( NCOMB, NSCI );
static_assert(NCOL >= CANVAS_LENGTH);

void beta_calc(std::string baseName = "", std::vector<int> ext = {}) {
	ROOT::EnableImplicitMT();

	std::vector<std::string> fileNames {};
	if(ext.size() == 0) 
		throw std::invalid_argument("Must supply atleast one file name.");
	
	fileNames.reserve(ext.size());
	for(int s : ext)
		fileNames.emplace_back(baseName + std::to_string(s) + ".root");

	const size_t NF = fileNames.size();
	int nf = 0;

	using Point  = std::pair<double,double>;
	using Points = std::vector<Point>;
	
	Points p{};
	
	std::vector<TCanvas*> c{};

	for(const auto& fileName : fileNames) {
		auto model = RNTupleModel::Create();
		auto frs = model->MakeField<RNFRSCal>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
	
		TH1P* h1_sci_e[NSCI];
		TH1P* h1_tof[NCOMB];

		for(int i=0; i<NSCI; ++i)
			h1_sci_e[i] = new TH1P(Form("M%d SCI%s QDC mean [QDC units]", nf, label[i]), ORGB{col[i]}, 600, 300, 4000);
		for(int i=0; i<NCOMB; ++i)
			h1_tof[i] = new TH1P(Form("M%d SCI%s-SCI%s ToF [ns]", nf, label[comb[i].first], label[comb[i].second]), ORGB{col[i]}, 8000, -100, 100);

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);
			
			for(int i=0; i<NSCI; ++i) {
				auto* h = h1_sci_e[i];
				const auto& sci = frs->sci[i];

				if(sci.hits.size() != 1) continue;
				(*h)->Fill(sci.E);
			}

			for(int i=0; i<NCOMB; ++i) {
				auto [i0,i1] = comb[i];
				const auto& sci0 = frs->sci[i0].hits;
				const auto& sci1 = frs->sci[i1].hits;
				
				if(sci0.size() != 1 or sci1.size() != 1) continue;
				
				double tof = sci1.back().t - sci0.back().t;
				(*h1_tof[i])->Fill(tof);
			}
		}
		c.emplace_back(new TCanvas(Form("c%d", nf), Form("c%d", nf), 1800, 1200));
		TCanvas* cb = c.back();
		cb->Divide(CANVAS_LENGTH, 2);
		for(int i=0; i<NSCI; ++i) {
			cb->cd(i+1);
			(*h1_sci_e[i])->Draw(); 
		}
		for(int i=0; i<NCOMB; ++i) {
			cb->cd(CANVAS_LENGTH + i+1);
			(*h1_tof[i])->Draw();
		}

		/* Objects owned back by gROOT. We lose direct handles. */
		++nf;
	}

}
