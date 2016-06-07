//
// Created by ayrat on yestoday.
//

#ifndef SDF_SYNTH_H
#define SDF_SYNTH_H


#include <string>
#include <unordered_map>


extern "C" {
#include <aiger.h>
#include <mtr.h>
#include <cudd.h>
};

#include <cuddObj.hh>
#include "myassert.hpp"
#include "Timer.hpp"


using namespace std;


//class Grapher;


class Synth {

public:
    Synth() {}
    // NOTE: time_limit_sec is used for heuristics (I won't stop on reaching it)
    bool run(const string &aiger_file_name, const string &output_file_name, unsigned time_limit_sec=3600);  // -> returns 'is realizable'
    ~Synth();


private:
    Synth(const Synth &other);
    Synth &operator=(const Synth &other);


private:
    Timer timer;
    unsigned time_limit_sec;
    Cudd cudd;
    aiger *aiger_spec;

    unordered_map<unsigned, BDD> transition_func;  // _aiger_unsigned_lit_ to bdd
    BDD init;
    BDD error;
    BDD win_region;
    BDD non_det_strategy;

    unordered_map<unsigned, BDD> bdd_by_aiger_unlit;   // this is specially for amba2match benchmarks

    unordered_map<unsigned, unsigned> cudd_by_aiger;   // mapping from aiger stripped literals to cudd indexes
    unordered_map<unsigned, unsigned> aiger_by_cudd;   // mapping from cudd indexes to aiger stripped literals

//    Grapher* grapher;   //pointer: let's keep class def secret from the compiler


private:
    BDD get_bdd_for_sign_lit(unsigned lit);

    vector<BDD> get_controllable_vars_bdds();
    vector<BDD> get_uncontrollable_vars_bdds();

    vector<BDD> get_bdd_vars(bool(*filter_func)(char *));

    void introduce_error_bdd();

    void compose_init_state_bdd();

    void compose_transition_vector();

    BDD pre_sys(BDD dst);  // also ensures that error is not violated

    BDD calc_win_region();

    BDD get_nondet_strategy();

    unordered_map<unsigned, BDD> extract_output_funcs();

    unsigned next_lit();

    unsigned get_optimized_and_lit(unsigned a_lit, unsigned b_lit);

    unsigned walk(DdNode *a_dd);

    void model_to_aiger(const BDD &c_signal, const BDD &func);

    vector<BDD> get_substitution();
};


#endif //SDF_SYNTH_H
