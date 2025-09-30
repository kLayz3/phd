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
			
			sci_raw.qdc[i] = static_cast<u16>( *sci_go4._qdc_raw[i] );
		}

		data->h1_sci_ml[s]->Fill( sci_raw.tdc[0].size() );	
		data->h1_sci_mr[s]->Fill( sci_raw.tdc[1].size() );	
	}

	for(int s=0; s < (int)data->tpc.size(); ++s) {
		auto& tpc_go4 = data->tpc[s];
		auto& tpc_raw = data->inner().tpc[s];

		/* Delay lines first. */
		for(int i=0; i<2; ++i) {
			Int_t nhits_l = *tpc_go4._tpc_ltn[i];
			Int_t nhits_r = *tpc_go4._tpc_ltn[i];
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
			auto& tdc_l = tpc_raw.tdc_l[i];
			auto& tdc_r = tpc_raw.tdc_r[i];
			tdc_l.resize(nhits_l);
			tdc_r.resize(nhits_r);
			std::memcpy(tdc_l.data(), tpc_go4._tpc_lt[i], nhits_l * sizeof(i32));
			std::memcpy(tdc_r.data(), tpc_go4._tpc_rt[i], nhits_r * sizeof(i32));	
		}

		/* Anodes (ADC + TDC) */
		for(int i=0; i<4; ++i) {
			Int_t nhits = *tpc_go4._tpc_atn[i];
			if(nhits < 0 || nhits > TFRSMapCont::TPC::MAX_SIZE) {
				WARN("Go4 reads: %d hits in TPC[%d] > %d or negative. Anode line[%d]. Dumping it.\n", 
					nhits, s, TFRSMapCont::TPC::MAX_SIZE, i);
				return;
			}
			auto& tdc = tpc_raw.tdc_a[i];
			tdc.resize(nhits);
			std::memcpy(tdc.data(), tpc_go4._tpc_at[i], nhits * sizeof(i32));
			tpc_raw.adc[i] = tpc_go4._tpc_aa[i]; 
		}

		Int_t nhits = *tpc_go4._sci_timerefn;
		auto& ref = tpc_raw.tdc_ref;
		ref.resize(nhits);
		std::memcpy(ref.data(), tpc_go4._sci_timeref, nhits * sizeof(i32));

		data->h1_tpc_ml[s] ->Fill( tpc_raw.tdc_l[0].size() );
		data->h1_tpc_mr[s] ->Fill( tpc_raw.tdc_r[1].size() );
		data->h1_tpc_ma1[s]->Fill( tpc_raw.tdc_a[1].size() );
		data->h1_tpc_ma2[s]->Fill( tpc_raw.tdc_a[2].size() );
		
		for(int a = 0; a<4; ++a) {
			if(int ndelay = a >> 1;
				tpc_raw.tdc_a[a].size() == 1 && 
				tpc_raw.tdc_r[ndelay].size() == 1 &&
				tpc_raw.tdc_l[ndelay].size() == 1)
			{
				data->h1_tpc_csum[s][a] -> Fill(
					tpc_raw.tdc_l[ndelay][0] + tpc_raw.tdc_r[ndelay][0] - 
					(tpc_raw.tdc_a[a][0] << 1)
				);
			}
		}
	}

	for(int s=0; s < (int)data->music.size(); ++s) {
		auto& music_go4 = data->music[s];
		auto& music_raw = data->inner().music[s];
	
		for(int i=0; i < music_raw.size(); ++i)
			music_raw.e[i] = static_cast<u16>(music_go4._music_raw[i]);
	}
}
