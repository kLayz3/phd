void extract_pad(
	const char* cname,
	int padn,
	TFile *f = gFile
) {
	if(!f) return;
	auto* c = dynamic_cast<TCanvas*>(f->Get(cname));
	if(!c) return;

	auto* pad = dynamic_cast<TPad*>(c->GetPad(padn));
	if(!pad) return;

	auto* out = new TCanvas("out", "Extracted pad", 1600, 1000);
	out->cd();

	auto* cpy = static_cast<TPad*>(pad->Clone("padc"));
	cpy->Draw();
	cpy->SetPad(0.,0., 1.,1.);
	cpy->Modified();

	out->Modified();
	out->Update();
}
