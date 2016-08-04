The safety synthesis tool submitted to the SYNTCOMP.
Dependencies:

- modifed aiger-1.9.4.
- cudd-3.0.0
- spdlog

To build, install dependencies into folder `<third_parties>`.
Modify `setup.sh` to provide the path to that folder.
Use the standard cmake routine for building:

- `. ./setup.sh`
- create build folder `<build>`
- `cd <build>`
- `cmake ..`
- `make`

Feel free to ask me questions, since the tool is raw.
