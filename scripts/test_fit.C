#include "../includes/PolyFitter.hxx"

double r();

void test_fit() {
	std::vector<double> x {-3,    -2,    -1,    0,    1,    2,    3};
	std::vector<double> y {-6.11, -4.20, -1.89, 0.11, 2.24, 3.89, 6.04};

	int N = (int)x.size();

	StaticPolyFitter<7,1> fit {x};

	auto [a0, a1] = fit.Fit(y);

	TGraph* g = new TGraph(N, x.data(), y.data());
	g->SetTitle("Static fit");
	g->SetMarkerStyle(20);   // e.g. full circle
	g->SetMarkerSize(2.5);   // default is ~1.0
							 //
	TF1* gFit = new TF1("f", "[0] + [1]*x", -4, 4);
	gFit->SetParameter(0, a0);
	gFit->SetParameter(1, a1);

	TCanvas* c = new TCanvas("c", "c", 1600, 1200);
	g->Draw("AP");
	gFit->Draw("SAME");


	/* Dynamic: */
	
	srand(time(NULL));
	N = (rand() % 20) + 6;
	printf("N-dyn: %d\n", N);
	x.clear(), y.clear();
	x.reserve(N), y.reserve(N);

	double a = 1.5 + (r()-0.5) / 5.0;
	double b = r() * 10;
	
	for(int i=0; i<N; ++i) {
		x.push_back(i);
		y.push_back( a*i + b + r()*2 );
	}
	auto [a0_d, a1_d] = PolyFit<1>(x, y);
	printf("Dfit: %.2f, %.2f\n", a0_d, a1_d);

	TGraph* gd = new TGraph(N, x.data(), y.data());
	gd->SetTitle("Dynamic fit");
	gd->SetMarkerStyle(20);   // e.gd. full circle
	gd->SetMarkerSize(2.5);   // default is ~1.0
							 //
	TF1* gFitD = new TF1("f", "[0] + [1]*x", x.front() - 1, x.back() + 1);
	gFitD->SetParameter(0, a0_d);
	gFitD->SetParameter(1, a1_d);

	TCanvas* cd = new TCanvas("cd", "cd", 1600, 1200);
	gd->Draw("AP");
	gFitD->Draw("SAME");
}

double r() {
	return static_cast<double>(rand()) / RAND_MAX;
}
