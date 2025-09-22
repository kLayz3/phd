#include "TFRSMapProc.h"
#include "libs.hh"
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
	data->Clear();

	for(int s=0; s < (int)data->sci.size(); ++s) {
		auto& sci = data->sci[s]; 
		
		for(int i=0; i<2; ++i) { // 0 => left, 1 => right
			Int_t nhits = *sci._nhit_raw[i];
			if(nhits < 0 || nhits > TFRSMapCont::Sci::MAX_SIZE) {
				WARN("Go4 reads: %d hits in Sci[%d] > %d or negative. Dumping it.\n", 
					nhits, s, TFRSMapCont::Sci::MAX_SIZE);
				return;
			}
			auto& vec = sci.tdc[i];

			vec.resize(nhits);
			std::memcpy(vec.data(),	sci._data_raw[i], nhits * sizeof(i32));
			
			sci.qdc[i] = static_cast<u16>( *sci._qdc_raw[i] );
		}

		data->h1_sci_ml[s]->Fill( sci.tdc[0].size() );	
		data->h1_sci_mr[s]->Fill( sci.tdc[1].size() );	
	}

	for(int s=0; s < (int)data->tpc.size(); ++s) {
		auto& tpc = data->tpc[s];

		/* Delay lines first. */
		for(int i=0; i<2; ++i) {
			Int_t nhits_l = *tpc._tpc_ltn[i];
			Int_t nhits_r = *tpc._tpc_ltn[i];
			if(nhits_l < 0 || nhits_l > TFRSMapCont::TPC::MAX_SIZE) {
				WARN("Go4 reads: %d hits in TPC[%d] > %d or negative. Left delay line[%d]. Dumping it.\n", 
					nhits_l, s, TFRSMapCont::TPC::MAX_SIZE, i);
				return;
			}
			if(nhits_r < 0 || nhits_r > TFRSMapCont::TPC::MAX_SIZE) {
				WARN("Go4 reads: %d hits in TPC[%d] > %d or negative. Right delay line[%d]. Dumping it.\n", 
					nhits_r, s, TFRSMapCont::TPC::MAX_SIZE, i);
				return;
			}
			auto& tdc_l = tpc.tdc_l[i];
			auto& tdc_r = tpc.tdc_r[i];
			tdc_l.resize(nhits_l);
			tdc_r.resize(nhits_r);
			std::memcpy(tdc_l.data(), tpc._tpc_lt[i], nhits_l * sizeof(i32));
			std::memcpy(tdc_r.data(), tpc._tpc_rt[i], nhits_r * sizeof(i32));	
		}

		/* Anodes (ADC + TDC) */
		for(int i=0; i<4; ++i) {
			Int_t nhits = *tpc._tpc_atn[i];
			if(nhits < 0 || nhits > TFRSMapCont::TPC::MAX_SIZE) {
				WARN("Go4 reads: %d hits in TPC[%d] > %d or negative. Anode line[%d]. Dumping it.\n", 
					nhits, s, TFRSMapCont::TPC::MAX_SIZE, i);
				return;
			}
			auto& tdc = tpc.tdc_a[i];
			tdc.resize(nhits);
			std::memcpy(tdc.data(), tpc._tpc_at[i], nhits * sizeof(i32));
			tpc.adc[i] = tpc._tpc_aa[i]; 
		}

		data->h1_tpc_ml[s] ->Fill( tpc.tdc_l[0].size() );
		data->h1_tpc_mr[s] ->Fill( tpc.tdc_r[1].size() );
		data->h1_tpc_ma1[s]->Fill( tpc.tdc_a[1].size() );
		data->h1_tpc_ma2[s]->Fill( tpc.tdc_a[2].size() );
	}

	for(auto& music : data->music) {
		for(int i=0; i < (int)LEN(music.e); ++i)
			music.e[i] = static_cast<u16>(music._music_raw[i]);
	}

	data->tpat = static_cast<u32>(*data->_pattern);
}
