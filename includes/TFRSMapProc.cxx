#include "TFRSMapProc.h"
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

	for(int i=0; i < (int)this->tpc.size(); ++i) {
		TFRSMapProc::TPC& tpc = this->tpc[i];
		tpc._tpc_aa = &sort->tpc_a[i][0];
		for(int j=0; j<2; ++j) {
			tpc._tpc_lt[j]  = &sort->tpc_lt[i][j][0];
			tpc._tpc_rt[j]  = &sort->tpc_rt[i][j][0];
			tpc._tpc_ltn[j] = &sort->tpc_nhit_lt[i][j];
			tpc._tpc_rtn[j] = &sort->tpc_nhit_rt[i][j];
		}
		for(int j=0; j<4; ++j) { 
			tpc._tpc_at[j]  = &sort->tpc_dt[i][j][0];
			tpc._tpc_atn[j] = &sort->tpc_nhit_dt[i][j];
		}
		tpc._sci_timerefn = &sort->tpc_nhit_timeref[i];
		tpc._sci_timeref = &sort->tpc_timeref[i][0];
	}

	this->music[0]._music_raw = &sort->music_e1[0];
	this->music[1]._music_raw = &sort->music_e2[0];
	this->_pattern = &sort->pattern;
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

	out.inner().tpat = static_cast<u32>(*_pattern);

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
		auto& tdc = tpc_raw.tdc; 

		/* Find how many elements we have to allocate. 
		 * This assignment is potentially slow. */

		for(int i=0; i<2; ++i) {
			_nhits_l[i] = *tpc_go4._tpc_ltn[i];
			_nhits_r[i] = *tpc_go4._tpc_rtn[i];
			if(_nhits_l[i] < 0 || _nhits_r[i] < 0) {
				WARN("Go4 reads: %d,%d < 0 (l,r) hits in TPC[%d], delay line[%d]. Skipping it. \n", 
					_nhits_l[i], _nhits_r[i], s, i);
				return;
			}
		}
		for(int i=0; i<4; ++i) {
			_nhits_a[i] = *tpc_go4._tpc_atn[i];
			if(_nhits_a[i] < 0) {
				WARN("Go4 reads: %d < 0 hits in TPC[%d], anode[%d].\n", _nhits_a[i], s, i);
				return;
			}
		}
		_nhits_s = *tpc_go4._sci_timerefn;
		if(_nhits_s < 0) {
			WARN("Go4 reads: %d < 0 hits in TPC[%d], sci ref.\n", _nhits_s, s);
			return;
		}
		if(s == 5 || s == 6) {
			if(_nhits_s > 0)
				WARN("Found a sci-ref hit in TPC%d - nhit = %d\n", s, _nhits_s);
		}
		
		/* We don't really care about more than RNTPCMap::MAX_SIZE elements. 
		 * Usually, that multiplicity indicates an inconsistent multi-hit event. */
		Int_t N_max = std::min (
			mnd::max (
				*std::max_element(std::begin(_nhits_l), std::end(_nhits_l)),
				*std::max_element(std::begin(_nhits_r), std::end(_nhits_r)),
				*std::max_element(std::begin(_nhits_a), std::end(_nhits_a)),
				_nhits_s
			), 
			RNTPCMap::MAX_SIZE
		);
		tdc.resize(N_max);
		std::fill_n(tdc.begin(), N_max, RNTPCMap::Measurement{});
		
		int _n_elems;
		
		/* Delay line left. */
		for(int i=0; i<2; ++i) {
			_n_elems = std::min(_nhits_l[i], N_max);
			for(int n=0; n < _n_elems; ++n) 
				_temp[n] = tpc_go4._tpc_lt[i][n];

			std::sort(_temp.begin(), _temp.begin() + _n_elems);

			for(int n=0; n < _n_elems; ++n) 
				tdc.at(n).tdc_l[i] = _temp[n];
		}
		/* Delay line right. */
		for(int i=0; i<2; ++i) {
			_n_elems = std::min(_nhits_r[i], N_max);
			for(int n=0; n < _n_elems; ++n) 
				_temp[n] = tpc_go4._tpc_rt[i][n];

			std::sort(_temp.begin(), _temp.begin() + _n_elems);

			for(int n=0; n < _n_elems; ++n) 
				tdc.at(n).tdc_r[i] = _temp[n];
		}
		/* Anodes (TDC + ADC) */
		for(int i=0; i<4; ++i) {
			_n_elems = std::min(_nhits_a[i], N_max);
			for(int n=0; n < _n_elems; ++n) 
				_temp[n] = tpc_go4._tpc_at[i][n];

			std::sort(_temp.begin(), _temp.begin() + _n_elems);

			for(int n=0; n < _n_elems; ++n) 
				tdc.at(n).tdc_a[i] = _temp[n];

			tpc_raw.adc[i] = tpc_go4._tpc_aa[i]; 
		}
		
		/* Sci reference. */
		_n_elems = std::min(_nhits_s, N_max);
		for(int n=0; n < _n_elems; ++n) 
			_temp[n] = tpc_go4._sci_timeref[n];

		std::sort(_temp.begin(), _temp.begin() + _n_elems);

		for(int n=0; n < _n_elems; ++n) 
			tdc.at(n).tdc_ref = _temp[n];

		/* ----------------------------------- */
		out.h1_tpc_ml[s] ->Fill( _nhits_l[0] );
		out.h1_tpc_mr[s] ->Fill( _nhits_r[1] );
		out.h1_tpc_ma1[s]->Fill( _nhits_a[1] );
		out.h1_tpc_ma2[s]->Fill( _nhits_a[2] );
	
		/* Only fill histos if N_max > 1. */
		if(N_max == 0) continue;

		for(int a = 0; a<4; ++a) { /* Loop over anode indices. */
			if(_nhits_a[a] != 1) continue;

			if(int ndelay = a >> 1;
				_nhits_r[ndelay] == 1 &&
				_nhits_l[ndelay] == 1)
			{
				out.h1_tpc_csum[s][a] -> Fill(
					tdc[0].tdc_l[ndelay] + tdc[0].tdc_r[ndelay] - 
					  2 * tdc[0].tdc_a[a]
				);
			}
			if(_nhits_s == 1)
			{
				out.h1_tpc_ydiff[s][a] -> Fill(
					tdc[0].tdc_a[a] -
					tdc[0].tdc_ref
				);
			}
		}
	} /* End loop over TPC's (s) */

	/* MUSIC's */
	for(int s=0; s < (int)music.size(); ++s) {
		const auto& music_go4 = this->music[s];
		auto& music_raw = out.inner().music[s];
	
		for(int i=0; i < music_raw.size(); ++i)
			music_raw.e[i] = static_cast<u16>(music_go4._music_raw[i]);
	}
}
