#include "PolyFitter.h"

template void PolyFit< 1>(const std::vector<double>& , const std::vector<double>& , std::array<double,  2>& );
template void PolyFit< 2>(const std::vector<double>& , const std::vector<double>& , std::array<double,  3>& );
template void PolyFit< 3>(const std::vector<double>& , const std::vector<double>& , std::array<double,  4>& );
template void PolyFit< 4>(const std::vector<double>& , const std::vector<double>& , std::array<double,  5>& );
template void PolyFit< 5>(const std::vector<double>& , const std::vector<double>& , std::array<double,  6>& );
template void PolyFit< 6>(const std::vector<double>& , const std::vector<double>& , std::array<double,  7>& );
template void PolyFit< 7>(const std::vector<double>& , const std::vector<double>& , std::array<double,  8>& );
template void PolyFit< 8>(const std::vector<double>& , const std::vector<double>& , std::array<double,  9>& );
template void PolyFit< 9>(const std::vector<double>& , const std::vector<double>& , std::array<double, 10>& );
template void PolyFit<10>(const std::vector<double>& , const std::vector<double>& , std::array<double, 11>& );
template void PolyFit<11>(const std::vector<double>& , const std::vector<double>& , std::array<double, 12>& );
template void PolyFit<12>(const std::vector<double>& , const std::vector<double>& , std::array<double, 13>& );

template std::array<double,  2> PolyFit< 1>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  3> PolyFit< 2>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  4> PolyFit< 3>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  5> PolyFit< 4>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  6> PolyFit< 5>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  7> PolyFit< 6>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  8> PolyFit< 7>(const std::vector<double>& , const std::vector<double>& );
template std::array<double,  9> PolyFit< 8>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 10> PolyFit< 9>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 11> PolyFit<10>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 12> PolyFit<11>(const std::vector<double>& , const std::vector<double>& );
template std::array<double, 13> PolyFit<12>(const std::vector<double>& , const std::vector<double>& );
