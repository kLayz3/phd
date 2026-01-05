#include "../includes/PolyFitter.hxx"

void test_fit(std::string fileName = "") {
	std::vector<double> x {-3,    -2,    -1,    0,    1,    2,    3};
	std::vector<double> y {-6.11, -4.20, -1.89, 0.11, 2.24, 3.89, 6.04};

	const int N =  (int)x.size();
	const int Ny = (int)y.size();
	assert(N == Ny);

	StaticPolyFitter<7,1> fit(x);

	auto [a0, a1] = fit.Fit(y);

	TGraph* g = new TGraph(N, x.data(), y.data());
	g->SetMarkerStyle(20);   // e.g. full circle
	g->SetMarkerSize(2.5);   // default is ~1.0
							 //
	TF1* gFit = new TF1("f", "[0] + [1]*x", -4, 4);
	gFit->SetParameter(0, a0);
	gFit->SetParameter(1, a1);

	TCanvas* c = new TCanvas("c", "c", 1600, 1200);
	g->Draw("AP");
	gFit->Draw("SAME");

	printf("A0: %.3f ; A1: %.3f\n", a0, a1);
}
