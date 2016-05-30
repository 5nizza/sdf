#include <fstream>
#include <iostream>

#include <spdlog/spdlog.h>

#include "Synth.hpp"


using namespace std;


void print_usage_exit(){
    cout << "sdf input_aiger_file [-o output_file]" << endl;
    cout << "-o output_file      is optional: if output_file=stdout," << endl
         << "                    then model is printed to stdout." << endl
         << "                    If not given, then model is calculated, but not output." << endl;
    exit(0);
}


int main (int argc, char *argv[]) {
    string output_file_name;
    string input_file_name;

    if (argc == 1 || !(argc == 4 || argc == 2))
        print_usage_exit();

    input_file_name = argv[1];

    if (argc == 4) {
        if (string(argv[2]) != "-o")
            print_usage_exit();

        output_file_name = argv[3];
    }

    // setup logging
    auto console = spdlog::stdout_logger_mt("console", true);
    spdlog::set_pattern("%H:%M:%S %v ");
    // end of logger setup

    Synth synthesizer;  //TODO: current: account for "" and "stdout"
    bool is_realizable = synthesizer.run(input_file_name, output_file_name);

    cout << (is_realizable ? "realizable":"unrealizable") << endl;
    return is_realizable? 10:20;
}
