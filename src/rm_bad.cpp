#include <fstream>
#include <iostream>

extern "C" {
#include <aiger.h>
};

#include "ArgsParser.h"


using namespace std;


void print_usage_exit() {
    cout << "I remove the 'bad' output and rename prefix controllable_ from outputs" << endl
         << "Usage:" << endl
         << "<tool> <input_aiger_file>" << endl;
    exit(0);
}


void main(int argc, char *argv[]) {
    ArgsParser parser(argc, argv, 2);
    if (parser.cmdOptionExists("-h") || argc < 2)
        print_usage_exit();

    string input_file_name = argv[1];
}