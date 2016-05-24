#include <fstream>
#include "Logger.hpp"
#include "Synth.hpp"



void print_usage_exit(){
    cout << "Run: <input_aiger_file> -o [output_file]" << endl;
    cout << "([output_file] is optional)" << endl;
    exit(0);
}


int main (int argc, char *argv[]) {
    if (argc == 1 || !(argc == 4 || argc == 2))
        print_usage_exit();

    if (argc == 4 && string(argv[2]) != "-o")
            print_usage_exit();

    string input_file_name = argv[1];
    string output_file_name = (argc == 4) ? argv[3]:"stdout";

    Synth synthesizer(input_file_name, output_file_name);
    bool is_realizable = synthesizer.run();

    return is_realizable? 10:20;
}
