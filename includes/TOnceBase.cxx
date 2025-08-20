#include "TOnceBase.h"

void TOnceBase::SetName(std::string name, const char* title)  { (void)title; this->_name = std::move(name); }
const char* TOnceBase::GetName() const { return this->_name.c_str(); }
