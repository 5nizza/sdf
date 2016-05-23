//
// Created by ayrat on today.
//

#ifndef SDF_BDDSYNTH_H
#define SDF_BDDSYNTH_H


extern "C" {
#include "aiger.h"
#include "cudd.h"
};

#include "cuddObj.hh"

#include <string>

using namespace std;


class BddSynth {

public:
    BddSynth(const string &aiger_file_name, const string &output_file_name);

    bool run(); // -> returns 'is realizable'


private:
    BddSynth(const BddSynth &other);

    BddSynth &operator=(const BddSynth &other);


private:
    Cudd *cudd;
    aiger *aiger_spec;
    aiger_symbol error_fake_latch;
    const string &_output_file_name;


public:
    virtual ~BddSynth();


private:
    BDD get_bdd_for_value(unsigned lit);

    vector<BDD> get_controllable_vars_bdds();

    vector<BDD> get_uncontrollable_output_bdds();

    vector<BDD> get_bdd_vars(int (*filter_func)(char *));

    int is_fake_error_lit(unsigned lit);

    void introduce_error_latch();

    BDD make_bdd_eq(BDD first, BDD second);

    BDD compose_init_state_bdd();

    BDD get_primed_variable_as_bdd(unsigned lit);

    BDD compose_transition_bdd();

    vector<BDD> get_all_latches_as_bdds();

    BDD prime_unprime_latches_in_bdd(BDD bdd, int should_prime);

    BDD prime_latches_in_bdd(BDD bdd);

    BDD unprime_latches_in_bdd(BDD bdd);

    BDD pre_sys_bdd(BDD dst_states, BDD transition);

    BDD calc_win_region(BDD init, BDD transition);

    BDD get_nondet_strategy(BDD win_region, BDD transition);

    vector<BDD> extract_output_funcs(BDD nondet_strategy);

    unsigned next_lit();

    unsigned get_optimized_and_lit(unsigned a_lit, unsigned b_lit);

    unsigned walk(DdNode *a_dd);

    void model_to_aiger(BDD &c_signal, BDD &func);
};


#endif //SDF_BDDSYNTH_H
