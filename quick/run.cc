#include "monad/monad.hxx"

#include <csignal>
#include <unistd.h>

#include "util/CLI.h"
#include "util/MacroHelpers.h"

int main(int argc, char* argv[]) {
	signal(SIGINT , mnd::sig_callback_handler);
	signal(SIGSEGV, mnd::sig_callback_handler);
	CLI::App app{"This program will delegate the arguments to the underlying program accordingly."};

	int nproc = 1;
	std::string prog, config, section="";
	bool help_underlying = false;
	
	add_logged_option(app, "-n,--nproc", nproc, "Number of subthreads " BOLD "(not implemented yet)" KNRM ".")
		->check(CLI::PositiveNumber)
		->each( [](const std::string& match) {
			int d = stoi(match);
			if( d & (d-1) ) {
				throw CLI::ValidationError("Number isn't a perfect power of 2.");
			}
		});
	add_logged_option<DisplayDefault::No>(app, "-p,--program", prog, "Program name")
		->required()
		->check(CLI::ExecPermissions);
	
	add_logged_option<DisplayDefault::No>(app, "-c,--config", config, "Config file name")
		->required()
		->check(CLI::ReadPermissions);

	add_logged_option(app, "-s,--section", section, "Section name. If left empty, the whole file becomes a single scope.");
	add_logged_flag(app, "-v", help_underlying, "Pass --help flag to the underlying process");

	CLI11_PARSE(app, argc, argv);
	
	/* Parse the config file. */
	const std::string cfg_contents = ParseFileToString(config);	
	std::string_view sec_contents;
	if(!section.empty()) {
		try {
			sec_contents = mnd::extract_section_body(cfg_contents, section).value();
		} catch(const std::exception& e) {
			ERROR("Contents empty, invalid or not section not found? Section name: \'%s\'\n", section.c_str());
		}
	} else {
		sec_contents = cfg_contents;
	}
	
	mnd::Argv fwd_argv = mnd::parse_argv(sec_contents, prog);
	if(help_underlying)
		fwd_argv.push_back("--help");

	WARN("Replacing the run script with \'%s\'. GLHF\n\n", fwd_argv[0]);
	execvp(fwd_argv[0], fwd_argv.data() );
	
	return 0;
}
