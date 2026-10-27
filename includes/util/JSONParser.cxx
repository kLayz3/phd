#include "JSONParser.h"
#include "FilesystemUtils.h"

#include <iostream>

using json = nlohmann::json;

json ParseJSON(const std::string& fileName) {
	auto text = mnd::fs::parse_file_to_string(fileName);
	try {
		return json::parse( text );
	}
	catch(const json::parse_error& e) {
		constexpr std::size_t context = 80;

		// e.byte is 1-based for parse errors.
		const std::size_t pos =
			e.byte > 0 ? e.byte - 1 : 0;

		const std::size_t begin =
			pos > context ? pos - context : 0;

		const std::size_t end =
			std::min(text.size(), pos + context);

		std::cerr
			<< "JSON parse error:\n"
			<< "  " << e.what() << '\n'
			<< "  byte: " << e.byte << '\n'
			<< "  context:\n"
			<< text.substr(begin, end - begin) << '\n';

		throw;
	}
}
