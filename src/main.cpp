#include <fstream>
#include "Logger.hpp"
#include "Synth.hpp"



void print_usage(){
    cout << "Run with two arguments: <input_aiger_file> [output_file]" << endl;
    cout << "([output_file] is optional)" << endl;
}


int main (int argc, char *argv[]) {
    if (argc <= 2) {
        print_usage();
        exit(0);
    }
    string input_file_name = argv[1];
    string output_file_name = (argc == 3) ? argv[2]:"stdout";

    Synth synthesizer(input_file_name, output_file_name);
    bool is_realizable = synthesizer.run();
    if (is_realizable)
        return 10;
    return 20;
}
