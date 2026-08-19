#pragma once

#include "monad/monad.hxx"
#include "TFRSCalCont.h"
#include "TFRSHitCont.h"

struct TFRSHitProc : TProcessor <
	TFRSHitCont
	(TFRSCalCont)
> {
	using Base = TProcessor<TFRSHitCont(TFRSCalCont)>;

	int s2_bt_tracking_mask;
	enum class S2BeforeTargetPositionV { all, any } s2_bt_posv;

	TFRSHitProc(TFRSHitCont& , const TFRSCalCont& , int );
	TFRSHitProc() = default;

	void ProcessEntry() noexcept;
	void ProcessS2BT() noexcept;
	void ProcessS2AT() noexcept;
    void ProcessS3() noexcept;

	/* Local buffer containers. */
	std::array<double, 2> fit_result;
	std::vector<double> x, y, zx, zy;
	
	/* EXPERT target parameters. */
	double z0, tar_width;
};
