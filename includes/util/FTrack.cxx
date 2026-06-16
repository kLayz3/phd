#include "FTrack.h"

std::ostream& operator<<(std::ostream& os, const FTrack& rhs) {
	return os << "L -> " << rhs.l << ", Q: " << rhs.q;
}

