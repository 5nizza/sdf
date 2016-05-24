#include "Synth.hpp"
#include "Logger.hpp"


#define IS_NEGATED(l) ((l & 1) == 1)
#define STRIP_LIT(lit) (lit & ~1)
#define NEGATED(lit) (lit ^ 1)


BDD Synth::get_bdd_for_value(unsigned lit) {
    /* lit is an AIGER variable index with a 'sign' */

    unsigned stripped_lit = STRIP_LIT(lit);
    BDD res;

    if (stripped_lit == 0) {
        res = cudd.bddZero();
    }
    else if (aiger_is_input(aiger_spec, stripped_lit) ||
             aiger_is_latch(aiger_spec, stripped_lit)) {
        res = cudd.bddVar(stripped_lit);  /*internal mapping of `interfaces'*/
    }
    else { // aiger_and
        aiger_and *and_ = aiger_is_and(aiger_spec, stripped_lit);
        res = get_bdd_for_value(and_->rhs0) & get_bdd_for_value(and_->rhs1);
    }

    return IS_NEGATED(lit) ? !res:res;
}


vector<BDD> Synth::get_bdd_vars(int(*filter_func)(char *)) {
    vector<BDD> result;
    for (unsigned i = 0; i < aiger_spec->num_inputs; i++) {
        aiger_symbol symbol = aiger_spec->inputs[i];
        if ((*filter_func)(symbol.name)) {
            BDD out_var_bdd = get_bdd_for_value(symbol.lit);
            result.push_back(out_var_bdd);
        }
    }

    return result;
}


int starts_with_controllable(char *name) { return 0 == string(name).find("controllable"); }

int not_starts_with_controllable(char *name) { return 0 != string(name).find("controllable"); }


vector<BDD> Synth::get_controllable_vars_bdds() {
    return get_bdd_vars(starts_with_controllable);
}


vector<BDD> Synth::get_uncontrollable_output_bdds() {
    return get_bdd_vars(not_starts_with_controllable);
}


void Synth::introduce_error_bdd() {
    error = get_bdd_for_value(aiger_spec->outputs[0].lit);
}


BDD Synth::make_bdd_eq(BDD first, BDD second) {
    return (first & second) | (~first & ~second);
}


void Synth::compose_init_state_bdd() { // Initial state is 'all latches are zero'
    L_INF("compose_init_state_bdd..");

    init = cudd.bddOne();

    for (unsigned i = 0; i < aiger_spec->num_latches; i++) {
        BDD latch_var = get_bdd_for_value(aiger_spec->latches[i].lit);
        init = init & ~latch_var;
    }
}


void Synth::compose_transition_vector() {
    L_INF("compose_transition_vector..");

    for (unsigned i = 0; i < aiger_spec->num_latches; ++i)
        transition_rel[aiger_spec->latches[i].lit] = get_bdd_for_value(aiger_spec->latches[i].next);
}


BDD Synth::pre_sys(BDD dst) {
    /**
    Calculate predecessor states of given states.

        ∀u ∃c ∃t': tau(t,u,c,t') & dst(t') & ~error(t,u,c)

    We use the direct substitution optimization (since t' <-> BDD(t,u,c)), thus:

        ∀u ∃c: (!error(t,u,c)  &  (dst(t)[t <- bdd_next_t(t,u,c)]))

    Note that we do not replace t variables in the error bdd.

    :return: BDD of the predecessor states
    **/

    dst = dst.VectorCompose(get_substitution());  // TODO: passing big vectors around (compress (i/2) indexing scheme?)

    BDD result = ~error & dst;  // TODO: use AndAbstract

    vector<BDD> controllable = get_controllable_vars_bdds();
    if (!controllable.empty()) {
        // ∃c ~error & dst[new_t]
        BDD vars_cube = cudd.bddComputeCube(controllable.data(), NULL, (int) controllable.size());   //TODO: cache this?
        result = result.ExistAbstract(vars_cube);
    }
    // else stays the same

    vector<BDD> uncontrollable = get_uncontrollable_output_bdds();
    if (!uncontrollable.empty()) {
        // ∀u ∃c ~error & dst[new_t]
        BDD uncontrollable_cube = cudd.bddComputeCube(uncontrollable.data(), NULL, (int) uncontrollable.size());  //TODO: cache this?
        result = result.UnivAbstract(uncontrollable_cube);
    }
    // else stays the same

    return result;
}

vector<BDD> Synth::get_substitution() {
    vector<BDD> substitution;   // TODO: cache this?

    for (unsigned i=0; i < (unsigned)cudd.ReadSize(); ++i)
        if (aiger_is_latch(aiger_spec, i))
            substitution.push_back(transition_rel[i]);
        else
            substitution.push_back(cudd.bddVar(i));

    return substitution;
}


BDD Synth::calc_win_region() {
    /** Calculate a winning region for the safety game: win = greatest_fix_point.X [not_error & pre_sys(X)]
        :return: BDD representing the winning region
    **/

    L_INF("calc_win_region..");

    BDD new_ = cudd.bddOne();
    while (1) {
        BDD curr = new_;

        new_ = pre_sys(curr);

        if (!(init <= new_))  // checking containment
            return cudd.bddZero();

        if (new_ == curr)
            return new_;
    }
}


BDD Synth::get_nondet_strategy(BDD win_region) {
    /**
    Get non-deterministic strategy from the winning region.
    If the system outputs controllable values that satisfy this non-deterministic strategy,
    then the system wins.
    I.e., a non-deterministic strategy describes for each state all possible plausible output values
    (below is assuming W excludes error states)

        strategy(t,u,c) = ∃t' W(t) & T(t,i,c,t') & W(t')

    But since t' <-> bdd(t,i,o), (and since we use error=error(t,u,c)), we use:

        strategy(t,u,c) = ~error(t,u,c) & W(t) & W(t)[t <- bdd_next_t(t,u,c)]

    :return: non-deterministic strategy bdd
    :note: The strategy is non-deterministic -- determinization is done later.
    **/

    L_INF("get_nondet_strategy..");

    return ~error & win_region & win_region.VectorCompose(get_substitution());
}



vector<BDD> Synth::extract_output_funcs(BDD nondet_strategy) {
    /** The result vector respects the order of the controllable variables **/

    L_INF("extract_output_funcs..");

    vector<BDD> models;

    vector<BDD> controls = get_controllable_vars_bdds();
    for (unsigned i = 0; i < controls.size(); ++i) {
        aiger_symbol *aiger_input = aiger_is_input(aiger_spec, STRIP_LIT(controls[i].NodeReadIndex()));
        L_INF("getting output function for " << aiger_input->name);

        // TODO: beautify: the current control is [0], others = [1..]
        swap(controls[0], controls[i]);
        BDD c = controls[0];

        BDD c_arena;
        if (controls.size()) {
            BDD cube = cudd.bddComputeCube(&controls[1], NULL, (int) (controls.size() - 1));
            c_arena = nondet_strategy.ExistAbstract(cube);
        }
        else { //special case of a single control
            c_arena = nondet_strategy;
        }
        // Now we have: c_arena(t,u,c) = ∃c_others: nondet(t,u,c)
        // (i.e., c_arena talks about this particular c, about t and u)

        BDD c_can_be_true = c_arena.Cofactor(c);
        BDD c_can_be_false = c_arena.Cofactor(~c);

        // We need to intersect with can_be_true to narrow the search.
        BDD c_must_be_true = ~c_can_be_false & c_can_be_true;
        BDD c_must_be_false = c_can_be_false & ~c_can_be_true;
        // Note that we cannot use `c_must_be_true = ~c_can_be_false`,
        // since the negation can cause including tuples (t,i,o) that violate non_det_strategy.

        // TODO: deref BDDs?

        BDD c_care_set = c_must_be_true | c_must_be_false;

        // We use 'restrict' operation, but we could also do just:
        // c_model = care_set -> must_be_true
        // ..but this is (probably) less efficient, since we cannot set c=1 if it is not in care_set, but we could.
        //
        // Restrict on the other side applies optimizations to find smaller bdd.
        // It cannot be expressed using boolean logic operations since we would need to say:
        // must_be_true = ite(care_set, must_be_true, "don't care")
        // and "don't care" cannot be expressed in boolean logic.

        // Restrict operation:
        //   on care_set: must_be_true.restrict(care_set) <-> must_be_true

        BDD c_model = c_must_be_true.Restrict(c_care_set);
        models.push_back(c_model);

        nondet_strategy = nondet_strategy & make_bdd_eq(c, c_model);

        // TODO: add optimization 'variables elimination'
    }

    swap(controls[0], controls[controls.size() - 1]);

    return models;
}


unsigned Synth::next_lit() {
    /* return: next possible to add to the spec literal */
    return (aiger_spec->maxvar + 1) * 2;
}


unsigned Synth::get_optimized_and_lit(unsigned a_lit, unsigned b_lit) {
    if (a_lit == 0 || b_lit == 0)
        return 0;

    if (a_lit == 1 && b_lit == 1)
        return 1;

    if (a_lit == 1)
        return b_lit;

    if (b_lit == 1)
        return a_lit;

    if (a_lit > 1 && b_lit > 1) {
        unsigned a_b_lit = next_lit();
        aiger_add_and(aiger_spec, a_b_lit, a_lit, b_lit);
        return a_b_lit;
    }

    MASSERT(0, "impossible");
}


/*
Walk given DdNode node (recursively).
If a given node requires intermediate AND gates for its representation, the function adds them.
    Literal representing given input node is `not` added to the spec.

:returns: literal representing input node
*/
unsigned Synth::walk(DdNode *a_dd) {  // TODO: add caching
    if (Cudd_IsConstant(a_dd))
        return (unsigned) (a_dd == cudd.bddOne().getNode());  // in aiger: 0 / 1 = False / True

    // get an index of variable,
    // all variables used in BDDs are also present in AIGER
    unsigned a_lit = Cudd_NodeReadIndex(a_dd);

    DdNode *t_bdd = Cudd_T(a_dd);
    DdNode *e_bdd = Cudd_E(a_dd);

    unsigned t_lit = walk(t_bdd);
    unsigned e_lit = walk(e_bdd);

    // ite(a_bdd, then_bdd, else_bdd)
    // = a*then + !a*else
    // = !(!(a*then) * !(!a*else))
    // -> in general case we need 3 more ANDs

    unsigned a_t_lit = get_optimized_and_lit(a_lit, t_lit);

    unsigned na_e_lit = get_optimized_and_lit(NEGATED(a_lit), e_lit);

    unsigned n_a_t_lit = NEGATED(a_t_lit);
    unsigned n_na_e_lit = NEGATED(na_e_lit);

    unsigned ite_lit = get_optimized_and_lit(n_a_t_lit, n_na_e_lit);

    unsigned res = NEGATED(ite_lit);
    if (Cudd_IsComplement(a_dd))
        res = NEGATED(res);

    return res;
}


/* Update aiger spec with a definition of `c_signal` */
void Synth::model_to_aiger(BDD &c_signal, BDD &func) {
    unsigned c_lit = c_signal.NodeReadIndex();

    unsigned func_as_aiger_lit = walk(func.getNode());

    aiger_redefine_input_as_and(aiger_spec, c_lit, func_as_aiger_lit, func_as_aiger_lit);
}


bool Synth::run() {
    L_INF("synthesize..");

    introduce_error_bdd();
    compose_init_state_bdd();
    compose_transition_vector();

//  time_t t;
//  struct tm * now;
//
//  t = time(0);   // get time now
//  now = localtime( & t );
//  cout << now->tm_hour << ":" << now->tm_min << endl;
//
//  std::cout << "REORDERING: ORIGINAL: SIZE: " << this->cudd.ReadNodeCount() << std::endl;
//
//  int count = 0;
//  while (count != this->cudd.ReadNodeCount()) {
//    count = this->cudd.ReadNodeCount();
//
//    this->cudd.ReduceHeap(CUDD_REORDER_SIFT_CONVERGE);
//    std::cout << "REORDERING: SIFT_CONV: SIZE: " << this->cudd.ReadNodeCount() << std::endl;
//
//    this->cudd.ReduceHeap(CUDD_REORDER_LINEAR_CONVERGE);
//    std::cout << "REORDERING: LIN_CONV: SIZE: " << this->cudd.ReadNodeCount() << std::endl;
//  }
//
//  // now add annealing
//  count = 0;
//  while (count != this->cudd.ReadNodeCount()) {
//    count = this->cudd.ReadNodeCount();
//
//    // very slow, but very good results!
//    this->cudd.ReduceHeap(CUDD_REORDER_ANNEALING);
//    std::cout << "REORDERING: ANNLNG: SIZE: " << this->cudd.ReadNodeCount() << std::endl;
//
//    this->cudd.ReduceHeap(CUDD_REORDER_SIFT_CONVERGE);
//    std::cout << "REORDERING: SIFT_CONV: SIZE: " << this->cudd.ReadNodeCount() << std::endl;
//
//    this->cudd.ReduceHeap(CUDD_REORDER_LINEAR_CONVERGE);
//    std::cout << "REORDERING: LIN_CONV: SIZE: " << this->cudd.ReadNodeCount() << std::endl;
//  }
//
////  this->cudd.ReduceHeap(CUDD_REORDER_EXACT);
////  std::cout << "after6: " << this->cudd.ReadNodeCount() << std::endl;
//
//  t = time(0);   // get time now
//  now = localtime( & t );
//  cout << now->tm_hour << ":" << now->tm_min << endl;
//
//  std::cout << "VARIABLE ORDER" << std::endl;
//  std::cout << this->cudd.OrderString() << std::endl;
//  exit(0);

    BDD win_region = calc_win_region();

    if (win_region.IsZero())
        return 0;

    BDD non_det_strategy = get_nondet_strategy(win_region);
    // win_region = cudd.bddZero();  // TODO: kill refs to non-needed BDDs
    vector<BDD> models = extract_output_funcs(non_det_strategy);
    vector<BDD> c_signals = get_controllable_vars_bdds();

    for (unsigned i = 0; i < models.size(); ++i)
        model_to_aiger(c_signals[i], models[i]);

    int res = (output_file_name == "stdout") ?  //TODO: magic constant
              aiger_write_to_file(aiger_spec, aiger_ascii_mode, stdout) :
              aiger_open_and_write_to_file(aiger_spec, output_file_name.c_str());
    MASSERT(res, "Could not write result file");

    return 1;
}


Synth::Synth(const string &aiger_file_name,
            const string &output_file_name) : output_file_name(output_file_name) {

    cudd.AutodynEnable(CUDD_REORDER_SIFT);

    aiger_spec = aiger_init();
    const char *err = aiger_open_and_read_from_file(aiger_spec, aiger_file_name.c_str());
    MASSERT(err == NULL, err);
}


Synth::~Synth() {
    aiger_reset(aiger_spec);
}
