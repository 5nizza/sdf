#include <fstream>
#include <iostream>

#include <spdlog/spdlog.h>

#include "Synth.hpp"


using namespace std;


void print_usage_exit(){
    cout << "Accepts one argument <input_aiger_file>. " << endl
         << "Thinks, then prints the result (UNREALIZABLE/REALIZABLE) followed by a model." << endl;
    exit(0);
}


int main (int argc, char *argv[]) {
    string input_file_name;

    if (argc != 2)
        print_usage_exit();

    input_file_name = argv[1];

//    if (argc == 4) {
//        if (string(argv[2]) != "-o")
//            print_usage_exit();
//
//        output_file_name = argv[3];
//    }

    // setup logging
    auto console = spdlog::stdout_logger_mt("console", false);
    spdlog::set_pattern("%H:%M:%S %v ");
    spdlog::set_level(spdlog::level::off);
    // end of logger setup

    Synth synthesizer;
    bool is_realizable = synthesizer.run(input_file_name, "stdout");   // TODO: stdout does not work..

    return is_realizable? 10:20;
}
