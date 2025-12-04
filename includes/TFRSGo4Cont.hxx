#pragma once

class TFRSSortEvent;

#include "TFRSSortEvent.h"
#include "monad/monad.hxx"

struct TFRSGo4Cont : TRawContainer<TFRSSortEvent> {
	TFRSGo4Cont() {
		this->SetName("FRSSortEvent.");
	}
};
