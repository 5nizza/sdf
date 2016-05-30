//
// Created by ayrat on today.
//

#ifndef SDF_SYNTH_H
#define SDF_SYNTH_H


#include <string>
#include <tr1/unordered_map>


extern "C" {
#include <aiger.h>
#include <mtr.h>
#include <cudd.h>
};

#include <cuddObj.hh>
#include "myassert.hpp"
#include "Timer.hpp"


using namespace std;


class Grapher;


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

    tr1::unordered_map<unsigned, BDD> transition_rel;  // _aiger_unsigned_lit_ to bdd
    BDD init;
    BDD error;

    Grapher* grapher;   //pointer: let's keep class def secret from the compiler


private:
    BDD get_bdd_for_value(unsigned lit);

    vector<BDD> get_controllable_vars_bdds();
    vector<BDD> get_uncontrollable_output_bdds();

    vector<BDD> get_bdd_vars(bool(*filter_func)(char *));

    void introduce_error_bdd();

    BDD make_bdd_eq(BDD first, BDD second);

    void compose_init_state_bdd();

    void compose_transition_vector();

    BDD pre_sys(BDD dst);  // also ensures that error is not violated

    BDD calc_win_region();

    BDD get_nondet_strategy(BDD win_region);

    vector<BDD> extract_output_funcs(BDD nondet_strategy);

    unsigned next_lit();

    unsigned get_optimized_and_lit(unsigned a_lit, unsigned b_lit);

    unsigned walk(DdNode *a_dd);

    void model_to_aiger(BDD &c_signal, BDD &func);

    vector<BDD> get_substitution();
};


#endif //SDF_SYNTH_H
