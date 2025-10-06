#pragma once

#include "Rtypes.h"

/**
 * Base class API for all types that will implement `void ProcessEntry()` and `Int_t Write()` methods.
 * It holds refs/ptrs to associated containers, or the containers should be movable.
 */

class TProcessor {
public:
	TProcessor() = default;
	TProcessor(const TProcessor& ) = delete;
	TProcessor& operator=(const TProcessor& other) = delete;

	TProcessor(TProcessor&& other) = default;
	TProcessor& operator=(TProcessor&& other) = default;

	virtual ~TProcessor() = default;

public:
	virtual Int_t Write();
};

static_assert(std::is_move_constructible<TProcessor>::value,  "Type <TProcessor> needs a move ctor.");
static_assert(std::is_move_assignable<TProcessor>::value, "Type <TProcessor> needs move assignment op.");
