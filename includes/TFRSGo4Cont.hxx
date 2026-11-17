#pragma once

#if defined(__GNUC__) && !defined(__clang__)
#	pragma GCC diagnostic push
#	pragma GCC diagnostic ignored "-Wunused-parameter"
#elif defined(__clang__)
#	pragma clang diagnostic push
#	pragma clang diagnostic ignored "-Wunused-parameter"
#endif

#include "TFRSSortEvent.h"

#if defined(__GNUC__) && !defined(__clang__)
#	pragma GCC diagnostic pop
#elif defined(__clang__)
#	pragma clang diagnostic pop
#endif

#include "monad/monad.hxx"

struct TFRSGo4Cont : TRawContainer<TFRSSortEvent> {
	TFRSGo4Cont() {
		this->SetName("FRSSortEvent.");
	}
};
