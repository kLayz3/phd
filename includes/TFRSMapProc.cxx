#include "TFRSMapProc.h"
#include "TFRSMapCont.h"
#include "libs.hh"
#include <algorithm>
#include <cstring>

void TFRSMapProc::ProcessEntry() noexcept {
	switch(do_analysis) {
		case 1: 
			_ProcessEntry();
		default:
			break;
	}
}

void TFRSMapProc::_ProcessEntry() noexcept {
	data->Clean();

	data->inner().tpat = static_cast<u32>(*data->_pattern);

	for(int s=0; s < (int)data->sci.size(); ++s) {
		auto& sci_go4 = data->sci[s];
		auto& sci_raw = data->inner().sci[s];
		
		for(int i=0; i<2; ++i) { // [0] => left, [1] => right
			Int_t nhits = *sci_go4._nhit_raw[i];
			if(nhits < 0 || nhits > TFRSMapCont::Sci::MAX_SIZE) {
				WARN("Go4 reads: %d hits in Sci[%d] > %d or negative. Dumping it.\n", 
					nhits, s, TFRSMapCont::Sci::MAX_SIZE);
				return;
			}
			auto& vec = sci_raw.tdc[i];

			vec.resize(nhits);
			std::memcpy(vec.data(),	sci_go4._data_raw[i], nhits * sizeof(i32));
			std::sort(vec.begin(), vec.end());
			sci_raw.qdc[i] = static_cast<u16>( *sci_go4._qdc_raw[i] );
		}

		data->h1_sci_ml[s]->Fill( sci_raw.tdc[0].size() );	
		data->h1_sci_mr[s]->Fill( sci_raw.tdc[1].size() );	
	}

	for(int s=0; s < (int)data->tpc.size(); ++s) {
		auto const& tpc_go4 = data->tpc[s];
		auto& tpc_raw = data->inner().tpc[s];
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
		
		/* We don't really care about more than RNTPCMap::MAX_SIZE elements. 
		 * Usually, that multiplicity indicates an inconsistent multi-hit event. */
		Int_t N_max = std::min (
			util::max (
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
		data->h1_tpc_ml[s] ->Fill( _nhits_l[0] );
		data->h1_tpc_mr[s] ->Fill( _nhits_r[1] );
		data->h1_tpc_ma1[s]->Fill( _nhits_a[1] );
		data->h1_tpc_ma2[s]->Fill( _nhits_a[2] );
	
		/* Only fill if N_max > 1. */
		if(N_max == 0) continue;

		for(int a = 0; a<4; ++a) { /* Loop over anode indices. */
			if(_nhits_a[a] != 1) continue;

			if(int ndelay = a >> 1;
				_nhits_r[ndelay] == 1 &&
				_nhits_l[ndelay] == 1)
			{
				data->h1_tpc_csum[s][a] -> Fill(
					tdc[0].tdc_l[ndelay] + tdc[0].tdc_r[ndelay] - 
					(tdc[0].tdc_a[a] << 1)
				);
			}
			if(_nhits_s == 1)
			{
				data->h1_tpc_ydiff[s][a] -> Fill(
					tdc[0].tdc_a[a] -
					tdc[0].tdc_ref
				);
			}
		}
	} /* End loop over TPC's (s) */

	for(int s=0; s < (int)data->music.size(); ++s) {
		auto& music_go4 = data->music[s];
		auto& music_raw = data->inner().music[s];
	
		for(int i=0; i < music_raw.size(); ++i)
			music_raw.e[i] = static_cast<u16>(music_go4._music_raw[i]);
	}
}
