/* Parse through a 12C files with corresponding brho's (S2-S3) to
 * calibrate S2 - S3 ToF. */

#include "ROOT/RNTupleModel.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "ROOT/RNTupleWriter.hxx"
#include "ROOT/RCanvas.hxx"
#include "ROOT/RNTupleDS.hxx"
#include "ROOT/RDataFrame.hxx"
#include "../../includes/util/JSONParser.h"
#include "../../includes/IonOptics.hxx"
#include "../../includes/util/PolyFitter.h"
#include "../../includes/util/PrettyHisto.hxx"
#include "../../includes/util/GaussFitMax.hxx"
#include "../../includes/util/FitDrawer.hxx"

using namespace ROOT;
using namespace ROOT::Experimental;
using json = nlohmann::json;

constexpr const char* label[] = {
	"21", "22", "23", "24"
}; constexpr int NSCI = static_cast<int>( sizeof(label)/sizeof(*label) );

#define SCI_22_I 1
#define SCI_31_I 2

using FileBrho = std::pair<std::string, double>;

void tof_cal (
	std::string param_file = "", 
	std::string common_part_prefix = "",
	std::vector<FileBrho> f = {},
	std::string common_part_suffix = "",
	std::array<double,3> dt_cut = {1000, -100, 100}
) {
	if(f.size() < 2)
		throw std::runtime_error("Must supply at least two files.");
	
	ROOT::EnableImplicitMT();
	
	json j = ParseJSON(param_file);
	json& j_s3 = j["FRS"]["S3"];
	
	FRSIdParam s3;
	UNROLL_JSON_PARAM(s3, j_s3, 5);

	std::vector<double> x, y;

	TCanvas* c = new TCanvas("dt", "Delta T", 1800, 1200);
	c->Divide(2,2);
	int i0 = 0;
	for(const auto& [file, brho] : f) {
		std::string fileName = common_part_prefix + file + common_part_suffix;

		auto model = RNTupleModel::Create();
		auto frs = model->MakeField<RNFRSCal>("FRS");
		auto ntuple = RNTupleReader::Open(std::move(model), "h103", fileName);
		const double beta = phy::Beta(6, 12, brho);
		const double x0 = 1.0 / beta;
		
		auto* h1_dt = new TH1P(Form("((h1_dt%d))Delta t [25 ps]@SCI31 - SCI22@TOF Point %s", i0, file.c_str()), kMagenta+i0, dt_cut[0], dt_cut[1], dt_cut[2]);

		for(auto entryId : *ntuple) {
			ntuple->LoadEntry(entryId);

			const auto& sci22 = frs->sci[SCI_22_I];
			const auto& sci31 = frs->sci[SCI_31_I];
			if(sci22.hits.size() != 1 or sci31.hits.size() != 1) continue;
			
			double dt = sci31.hits[0].t - sci22.hits[0].t;
			h1_dt->Fill(dt);
		}

		auto [res, err_] = GaussFitMax(*h1_dt);
		x.push_back(x0);
		y.push_back(res[1]);

		c->cd(++i0);
		h1_dt->Draw();
	}
	std::array<double, 2> r;	
	auto graph = FitAndDraw<1>(x, y, r);
	
	TCanvas* cf = new TCanvas("cf", "Fit", 800, 500);
	graph.first->Draw("AP");
	graph.second->Draw("L SAME");
	gPad->SetGrid();

	std::cout << "x: " << x <<  ", y: " << y << std::endl;
	std::cout << setprecision(10) << "Final result: " << r[0] << ", " << r[1] << std::endl;
	
}
