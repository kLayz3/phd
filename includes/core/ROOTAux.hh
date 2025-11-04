#pragma once

#include <system_error>
#include <unistd.h>
#include <filesystem>
#include <cstdlib>
#include "libs.hh"

/* Whoever is including this, will probably talk to also `gSystem` handle. */
#include "TSystem.h" 

namespace util {

inline void CreateClingCache() {
	std::string cache_dir = 
		std::string(getenv("HOME") ? getenv("HOME") : "/tmp")
		+ "/.cling_cache_"
		+ std::string(getenv("HOSTNAME") ? getenv("HOSTNAME") : "11h_Unknown")
		+ "_p"
		+ std::to_string(getpid());

	std::error_code ec;
	setenv("CLING_MODULE_CACHE_PATH", cache_dir.c_str(), 1);
	std::filesystem::create_directories(cache_dir, ec);
	
	if(ec) {
		ERROR("Exception: %s . Failed to create cache dir at: " EMPH(%s), ec.message().c_str(), cache_dir.c_str());
	}

	if(!std::filesystem::is_directory(cache_dir)) {
		ERROR("Cache dir: %s created somehow, but it is not a directory?", cache_dir.c_str());
	}
	
	setenv("ROOT_DISABLE_MODULES","1",1);
	WARN("Successfully created cache dir: %s\n", cache_dir.c_str());
}

}
