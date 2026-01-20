/* There are certain forward declarations in MONAD that maybe users wish to
 * declare before the definitions passed along from MONAD itself.
 * 
 * An example would be the overloaded <<operator which recurses into array<?,?>|vector<?>
 * The non-container type `T` then needs to have `ostream& operator<<(ostream&, const T&)` symbol
 * visible before the general definition for array/vector overloads. 
 * So what? Just forward declare `struct T` and `ostream& operator<<(ostream&, const T&)` symbols?
 * Yes - but what if `T` is a nested type. Can't forward declare then. 
 *
 * Decoupling the declarations after general MONAD part I don't want to do since we might use it
 * for debugging arrays/vectors inside MONAD itself. 
 *
 * In that case, just include this forwarding file first, then define the (nested) structure, and include `monad.hxx`
 * after. */

/* Careful with this one, a lot of templates are in global namespace. */
#pragma once

#include <iostream>
#include <vector>
#include <array>

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& );
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& );
/* ^^^ Fwd declared for symbol visiblity in the function below
 * All underlying 'bare' `T` must have the overloaded operator defined at this point. */

template<typename T>
std::ostream& mnd_output_homogeneous_range_(std::ostream& os, const T* p, const std::size_t N) {
	os << '[';
	if(N > 0) os << p[0];
	for(std::size_t i = 1; i < N; ++i) {
		os << ", " << p[i];
	}
	os << ']';
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	return mnd_output_homogeneous_range_(os, v.data(), v.size());
}
template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& v) {
	return mnd_output_homogeneous_range_(os, v.data(), v.size());
}
