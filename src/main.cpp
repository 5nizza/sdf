#include <fstream>
#include <iostream>

#include <spdlog/spdlog.h>

#include "Synth.hpp"
#include "ArgsParser.h"


using namespace std;


void print_usage_exit() {
    cout << "<tool> <input_aiger_file> [-f]" << endl;
    cout << "Use flag -f if you need pure full models (with outputs defined and `bad` output removed)" << endl;
    exit(0);
}


int main (int argc, char *argv[]) {
    ArgsParser parser(argc, argv, 2);
    if (parser.cmdOptionExists("-h") || argc < 2)
        print_usage_exit();

    string input_file_name = argv[1];

    bool print_full_model = parser.cmdOptionExists("-f");

    string output_file_name = parser.cmdOptionExists("-o") ? parser.getCmdOption("-o") : "stdout";

    // setup logging
    auto console = spdlog::stdout_logger_mt("console", false);
    spdlog::set_pattern("%H:%M:%S %v ");
    //spdlog::set_level(spdlog::level::off);
    // end of logger setup

    Synth synthesizer(input_file_name, output_file_name, print_full_model, 3600);  // TODO: stdout does not work..
    bool is_realizable = synthesizer.run();

    return is_realizable? 10:20;
}
