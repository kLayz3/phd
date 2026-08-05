#include "FitSpline.h"

/* Explicitly instantiate the runtime polynomial degree routines. */

template std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, TGraphErrors*> 
mnd::detail::FitSplineImpl<fit_info::PROFILE_MAX> (TH2D* , double& , double& , double , uint32_t, Verbosity );

template std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, TGraphErrors*> 
mnd::detail::FitSplineImpl<fit_info::GAUSS_MAX> (TH2D* , double& , double& , double , uint32_t, Verbosity );

template std::tuple<std::vector<double>, TGraphErrors*, TGraph*> 
FitSpline<fit_info::PROFILE_MAX>(std::size_t, TH2D* , double , double , int , double , uint32_t, Verbosity );

template std::tuple<std::vector<double>, TGraphErrors*, TGraph*> 
FitSpline<fit_info::GAUSS_MAX>(std::size_t, TH2D* , double , double , int , double , uint32_t, Verbosity );
