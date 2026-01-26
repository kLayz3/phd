#include "TFRSMapProc.h"
#include "TFRSMapCont.h"
#include <algorithm>

#if defined(__GNUC__)
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Winvalid-offsetof"
#	pragma GCC diagnostic pop
#endif

void TFRSMapProc::SetupPointers() {
	TFRSGo4Cont & input = std::get<0>( this->in );
	if(input.raw() == nullptr) ERROR("TFRSGo4Cont, underlying pointer is null?");

	TFRSSortEvent* sort = input.raw();
#define MAP_SCI(x, SCI_LABEL) \
	if(x < N_VALID_SCI) { \
		this->sci[x]._nhit_raw[0] = &sort->tdc_nhit_sc##SCI_LABEL##l; \
		this->sci[x]._nhit_raw[1] = &sort->tdc_nhit_sc##SCI_LABEL##r; \
		this->sci[x]._data_raw[0] = &sort->tdc_sc##SCI_LABEL##l[0]; \
		this->sci[x]._data_raw[1] = &sort->tdc_sc##SCI_LABEL##r[0]; \
		this->sci[x]._qdc_raw[0]  = &sort->de_##SCI_LABEL##l; \
		this->sci[x]._qdc_raw[1]  = &sort->de_##SCI_LABEL##r; \
	}
	
	MAP_SCI(0, 21);
	MAP_SCI(1, 22);
	MAP_SCI(2, 31);
	MAP_SCI(3, 41);

	for(int s=0; s < (int)this->tpc.size(); ++s) {
		TFRSMapProc::TPC& tpc = this->tpc[s];
		tpc._tpc_aa = &sort->tpc_a[s][0];
		for(int i : {0,1}) { // Delay line index.
			tpc._tpc_ltn[i] = &sort->tpc_nhit_lt[s][i];
			tpc._tpc_rtn[i] = &sort->tpc_nhit_rt[s][i];
			tpc._tpc_lt [i] = &sort->tpc_lt[s][i][0];
			tpc._tpc_rt [i] = &sort->tpc_rt[s][i][0];
			
			for(int a : {0,1}) { // Anode index, underlying the delay line `i`.
				tpc._tpc_atn[i][a] = &sort->tpc_nhit_dt[s][2*i + a];
				tpc._tpc_at [i][a] = &sort->tpc_dt[s][2*i + a][0];
			}
		}
		tpc._sci_timerefn = &sort->tpc_nhit_timeref[s];
		tpc._sci_timeref = &sort->tpc_timeref[s][0];
	}

	this->music[0]._music_raw = &sort->music_e1[0];
	this->music[1]._music_raw = &sort->music_e2[0];
}


void TFRSMapProc::ProcessEntry() noexcept {
	switch(do_analysis) {
		case DoAnalysis::YES: 
			_ProcessEntry();
		case DoAnalysis::NO:
			break;
	}
}

void TFRSMapProc::_ProcessEntry() noexcept {
	SetupPointers();
	out.Clean();

	/* Scintillators. */
	for(int s=0; s < (int)sci.size(); ++s) {
		const auto& sci_go4 = this->sci[s];
		auto& sci_raw = out.inner().sci[s];
		auto& tdc = sci_raw.tdc;
	
		int nhits[2] = {0};      // [0] => left, [1] => right
		for(int i=0; i<2; ++i) { // [0] => left, [1] => right
			nhits[i] = *sci_go4._nhit_raw[i];
			if(nhits[i] < 0 || nhits[i]> RNSciMap::MAX_SIZE) {
				WARN("Go4 reads: %d hits in Sci[%d] > %d or negative. Dumping it.\n", 
					nhits[i], s, RNSciMap::MAX_SIZE);
				return;
			}
			sci_raw.qdc[i] = static_cast<u16>( *sci_go4._qdc_raw[i] );
		}
		Int_t N_max = std::min (
			std::max(nhits[0], nhits[1]),
			RNSciMap::MAX_SIZE
		);

		tdc.resize(N_max);
		std::fill_n(tdc.begin(), N_max, RNSciMap::Measurement{});
	
		std::array<Int_t, RNSciMap::MAX_SIZE> _tmp;
		int _n_elems;

		for(int i=0; i<2; ++i) { // [0] => left, [1] => right
			_n_elems = std::min(nhits[i], N_max);
			for(int n=0; n < _n_elems; ++n) 
				_tmp[n] = sci_go4._data_raw[i][n];
			
			std::sort(_tmp.begin(), _tmp.begin() + _n_elems);
			
			if(i == 0) for(int n=0; n < _n_elems; ++n) tdc[n].tdc_l = _tmp[n];
			else       for(int n=0; n < _n_elems; ++n) tdc[n].tdc_r = _tmp[n];
		}

		out.h1_sci_ml[s]->Fill( nhits[0] );	
		out.h1_sci_mr[s]->Fill( nhits[1] );	
		if(nhits[0] == nhits[1])
			out.h1_sci_diff_lr[s]->Fill( tdc[0].tdc_l - tdc[0].tdc_r );
	}

	/* TPC's */
	for(int s=0; s < (int)tpc.size(); ++s) {
		const auto& tpc_go4 = this->tpc[s];
		auto& tpc_raw = out.inner().tpc[s];
		auto& tdc_ref = tpc_raw.tdc_ref; // vec<i32>
		auto& tdc_dl  = tpc_raw.tdc;     // array<vec<Measurement>> 

		/* Find how many elements we have to allocate. 
		 * This assignment is potentially slow. */
		
		int _n_elems;

		_nhits_s = *tpc_go4._sci_timerefn;
		if(_nhits_s < 0) {
			WARN("Go4 reads: %d < 0 hits in TPC[%d], sci ref.\n", _nhits_s, s);
			return;
		}
		
		/* Sci reference. */
		_n_elems = std::min(_nhits_s, RNTPCMap::MAX_SIZE);
		tdc_ref.resize(_n_elems);
		for(int n=0; n < _n_elems; ++n) 
			_temp[n] = tpc_go4._sci_timeref[n];

		std::sort(_temp.begin(), _temp.begin() + _n_elems);

		for(int n=0; n < _n_elems; ++n) 
			tdc_ref.at(n) = _temp[n];
	
		/* Loop over delay-lines, sort out the data. */
		for(int i : {0,1}) {
			std::vector<RNTPCMap::Measurement>& tdc = tdc_dl[i]; 

			_nhits_l[i] = *tpc_go4._tpc_ltn[i];
			_nhits_r[i] = *tpc_go4._tpc_rtn[i];
			if(_nhits_l[i] < 0 || _nhits_r[i] < 0) {
				WARN("Go4 reads: %d,%d < 0 (l,r) hits in TPC[%d], delay line[%d]. Skipping it. \n", 
					_nhits_l[i], _nhits_r[i], s, i);
				return;
			}
			for(int a : {0,1}) {
				_nhits_a[i][a] = *tpc_go4._tpc_atn[i][a];
				if(_nhits_a[i][a] < 0) {
					WARN("Go4 reads: %d < 0 hits in TPC[%d], anode[%d].\n", _nhits_a[i][a], s, 2*i+a);
					return;
				}
			}

			/* We don't really care about more than RNTPCMap::MAX_SIZE elements. 
			 * Usually, that multiplicity indicates an inconsistent multi-hit event. */
			Int_t N_max = std::min (
				mnd::max (
					_nhits_l[i], _nhits_r[i],
					_nhits_a[i][0], _nhits_a[i][1]
				), 
				RNTPCMap::MAX_SIZE
			);
			tdc.resize(N_max);
			std::fill_n(tdc.begin(), N_max, RNTPCMap::Measurement{});
			
			/* Delay line left. */
			_n_elems = std::min(_nhits_l[i], N_max);
			for(int n=0; n < _n_elems; ++n) 
				_temp[n] = tpc_go4._tpc_lt[i][n];

			std::sort(_temp.begin(), _temp.begin() + _n_elems);

			for(int n=0; n < _n_elems; ++n)
				tdc[n].tdc_l = _temp[n];
			/* =================== */

			/* Delay line right. */
			_n_elems = std::min(_nhits_r[i], N_max);
			for(int n=0; n < _n_elems; ++n) 
				_temp[n] = tpc_go4._tpc_rt[i][n];

			std::sort(_temp.begin(), _temp.begin() + _n_elems);

			for(int n=0; n < _n_elems; ++n) 
				tdc[n].tdc_r = _temp[n];
			/* =================== */

			/* Anodes (TDC) */
			for(int a : {0,1}) {
				_n_elems = std::min(_nhits_a[i][a], N_max);
				for(int n=0; n < _n_elems; ++n)
					_temp[n] = tpc_go4._tpc_at[i][a][n];

				std::sort(_temp.begin(), _temp.begin() + _n_elems);

				for(int n=0; n < _n_elems; ++n)
					tdc[n].tdc_a[a] = _temp[n];
			}
			/* =================== */

			if(i == 0) {
				out.h1_tpc_ml[s] ->Fill( _nhits_l[0] );
				out.h1_tpc_mr[s] ->Fill( _nhits_r[1] );
				out.h1_tpc_ma1[s]->Fill( _nhits_a[0][0] );
				out.h1_tpc_ma2[s]->Fill( _nhits_a[0][1] );
			}
		
			if(_nhits_a[i][0] == 1 and _nhits_a[i][1] == 1) {
				out.h1_tpc_adiff[s][i] -> Fill(
					tdc.at(0).tdc_a[0] - tdc.at(0).tdc_a[1]
				);
			}

			if(_nhits_l[i] != 1 || _nhits_r[i] != 1) continue;
			
			for(int a : {0,1} ) {
				if(_nhits_a[i][a] == 1) {
					out.h1_tpc_csum[s][2*i + a] -> Fill(
						tdc[0].tdc_l + tdc.at(0).tdc_r - 
						2 * tdc[0].tdc_a[a]
					);

					out.h1_tpc_araw[s][2*i+a]->Fill(tdc.at(0).tdc_a[a]);

					if(_nhits_s == 1)
						out.h1_tpc_ydiff[s][2*i + a] -> Fill(
							tdc[0].tdc_a[a] -
							tdc_ref[0]
						);
				}
			}
		} // end of loop over delay-lines {0,1}
		
		if(_nhits_l[0] == 1 && _nhits_l[1] == 1) {
			out.h1_tpc_ldiff[s]->Fill (
				tdc_dl[0].at(0).tdc_l - tdc_dl[1].at(0).tdc_l
			);
			for(int d: {0,1} )
				out.h1_tpc_dl_lraw[s][d]->Fill(tdc_dl[d].at(0).tdc_l);
		}
		if(_nhits_r[0] == 1 && _nhits_r[1] == 1) {
			out.h1_tpc_rdiff[s]->Fill (
				tdc_dl[0].at(0).tdc_r - tdc_dl[1].at(0).tdc_r
			);
			for(int d: {0,1} )
				out.h1_tpc_dl_rraw[s][d]->Fill(tdc_dl[d].at(0).tdc_r);
		}
		
		/* Anodes (ADC) */
		for(int i=0; i<4; ++i)
			tpc_raw.adc[i] = tpc_go4._tpc_aa[i]; 

	} /* End loop over TPC's (s) */

	/* MUSIC's */
	for(int s=0; s < (int)music.size(); ++s) {
		const auto& music_go4 = this->music[s];
		auto& music_raw = out.inner().music[s];
	
		for(int i=0; i < music_raw.size(); ++i)
			music_raw.e[i] = static_cast<u16>(music_go4._music_raw[i]);
	}
}

void TTrigMapProc::ProcessEntry() noexcept {
	if(do_analysis == TFRSMapProc::DoAnalysis::NO) return;
		
	TFRSGo4Cont& input = std::get<0>( this->in );
	TFRSSortEvent* sort = input.raw();
	if(sort == nullptr) ERROR("TFRSGo4Cont, underlying pointer is null?");
	
	this->out.Clean();
	RNTrigMap& out = this->out.inner(); 
	
	u64 wr = sort->frs_wr; 
	if(wr == 0) return;	
	
	/* First assign tpat, because that's the nullable indicator. */
	u16 tpat = sort->pattern & 0xffff;
	out.tpat = tpat; 
	this->out.h1_full_tpat->Fill(tpat);
	
	out.wr.assign( static_cast<u32>(wr) );
	this->out.h1_wr_diff->Fill( out.DeltaT().value_or(NAN) / 1000. );

	/* Butcher the (local) tpat now, look at all the bits. */
	if(tpat == 0) {
		this->out.h1_tpat->Fill(0);
	} else {
		for(int i=1; i<=16; ++i) {
			if(tpat & 1)
				this->out.h1_tpat->Fill(i);
			tpat >>= 1;
		}
	}
}
