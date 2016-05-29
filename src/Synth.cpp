#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <spdlog/spdlog.h>
#include <sstream>
#include <vector>

#include "Synth.hpp"


#define IS_NEGATED(l) ((l & 1) == 1)
#define STRIP_LIT(lit) (lit & ~1)
#define NEGATED(lit) (lit ^ 1)


#define L_INF(message) {spdlog::get("console")->info() << message;}


class Grapher {
public:
    Grapher(aiger* spec) {
        this->spec = spec;
    }

    aiger* spec;
    unordered_map<unsigned, unordered_set<unsigned> > deps;
    unordered_map<unsigned, unsigned> freq_map;

    unordered_set<unsigned>
    get_deps(unsigned lit, aiger* spec) {
        auto stripped_lit = STRIP_LIT(lit);

        if (deps.find(stripped_lit) == deps.end()) {
            auto and_ = aiger_is_and(spec, stripped_lit);                       MASSERT(and_, "impossible");
            auto rhs0_deps = get_deps(and_->rhs0, spec);
            auto rhs1_deps = get_deps(and_->rhs1, spec);

            unordered_set<unsigned> a_deps;
            a_deps.insert(rhs0_deps.begin(), rhs0_deps.end());
            a_deps.insert(rhs1_deps.begin(), rhs1_deps.end());

            deps[stripped_lit] = a_deps;
        }

        return deps[stripped_lit];
    }

    void compute_deps() {
        deps[0] = unordered_set<unsigned>();  // value 'false' (need to insert this value, since we later call 'find')

        for (unsigned i=0; i < spec->num_latches; ++i)
            deps[spec->latches[i].lit].insert(spec->latches[i].lit);

        for (unsigned i=0; i < spec->num_inputs; ++i)
            deps[spec->inputs[i].lit].insert(spec->inputs[i].lit);

        // above we initialized values for latches and inputs ('real' variables)
        // now we calculate dependencies between them

        for (unsigned i=0; i < spec->num_ands; ++i)
            deps[spec->ands[i].lhs] = get_deps(spec->ands[i].lhs, spec);

        // at this point deps of all 'real' vars are computed:
        // - input vars do not have deps,
        // - latches vars are computed because dependencies of their next literals are computed
        compute_freq_map();
    }

    void compute_freq_map() {
        for (unsigned i=0; i < spec->num_latches; ++i) {
            auto latch_deps = deps[STRIP_LIT(spec->latches[i].next)];
            for (auto const & elem : latch_deps)
                ++freq_map[elem];
        }

        auto output_deps = deps[STRIP_LIT(spec->outputs[0].lit)];
        for (auto const & elem : output_deps)
            ++freq_map[elem];
    }

    // TODO:
//    unsigned is_referenced_only_once(unsigned and_lit) { }
//
//     TODO:
//    void distance(unsigned lit1, unsigned lit2) { }
//
    void dump_dot() {
        cout << "digraph latch_graph { rankdir=BT;" << endl;

        for (unsigned i=0; i < spec->num_inputs; i++)
            cout << spec->inputs[i].lit << "[shape=triangle];" << endl;
        for (unsigned i=0; i < spec->num_latches; i++)
            cout << spec->latches[i].lit << "[shape=diamond, color=magenta];" << endl;
        cout << STRIP_LIT(spec->outputs[0].lit) << "[shape=triangle, color=blue];" << endl;

        for (unsigned i=0; i < spec->num_latches; i++)
            for (auto const & src : deps[STRIP_LIT(spec->latches[i].next)])
                cout << src << "->" << spec->latches[i].lit << ";" << endl;

        for (auto const & src : deps[STRIP_LIT(spec->outputs[0].lit)])
            cout << src << "->" << STRIP_LIT(spec->outputs[0].lit) << ";" << endl;

        cout << "}" << endl;
    }
//    void dump_order(string file_name) { }
};



BDD Synth::get_bdd_for_value(unsigned lit) {
    /* lit is an AIGER variable index with a 'sign' */

    unsigned stripped_lit = STRIP_LIT(lit);
    BDD res;

    if (stripped_lit == 0) {
        res = cudd.bddZero();
    }
    else if (aiger_is_input(aiger_spec, stripped_lit) ||
             aiger_is_latch(aiger_spec, stripped_lit)) {
        res = cudd.ReadVars(stripped_lit/2);  // internal cudd mapping (we don't need to create the var, so I use ReadVars)
        MASSERT(res.NodeReadIndex() == stripped_lit/2, "that bug again: impossible: " << res.NodeReadIndex() << " vs " << stripped_lit/2 );
    }
    else { // aiger_and
        aiger_and *and_ = aiger_is_and(aiger_spec, stripped_lit);
        res = get_bdd_for_value(and_->rhs0) & get_bdd_for_value(and_->rhs1);
    }

    return IS_NEGATED(lit) ? (~res):res;
}


vector<BDD> Synth::get_bdd_vars(bool(*filter_func)(char *)) {
    vector<BDD> result;
    for (unsigned i = 0; i < aiger_spec->num_inputs; ++i) {
        aiger_symbol symbol = aiger_spec->inputs[i];
        if ((*filter_func)(symbol.name)) {
            BDD out_var_bdd = get_bdd_for_value(symbol.lit);
            result.push_back(out_var_bdd);
        }
    }

    return result;
}


bool starts_with_controllable(char *name) { return 0 == string(name).find("controllable"); }

bool not_starts_with_controllable(char *name) { return string::npos == string(name).find("controllable"); }


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
    return ~(first.Xor(second));
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

    dst = dst.VectorCompose(get_substitution());

    vector<BDD> controllable = get_controllable_vars_bdds();

    MASSERT(!controllable.empty(), "");

    // ∃c ~error & dst[new_t]
    BDD vars_cube = cudd.bddComputeCube(controllable.data(), NULL, (int) controllable.size());
    BDD result = dst.AndAbstract(~error, vars_cube);

    vector<BDD> uncontrollable = get_uncontrollable_output_bdds();
    if (!uncontrollable.empty()) {
        // ∀u ∃c ~error & dst[new_t]
        BDD uncontrollable_cube = cudd.bddComputeCube(uncontrollable.data(), NULL, (int) uncontrollable.size());
        result = result.UnivAbstract(uncontrollable_cube);
    }

    return result;
}

vector<BDD> Synth::get_substitution() {
    vector<BDD> substitution;

    substitution.push_back(cudd.ReadVars(0));   // special variable
    for (unsigned i=1; i < (unsigned)cudd.ReadSize(); ++i) {
        if (aiger_is_latch(aiger_spec, i*2))
            substitution.push_back(transition_rel.find(i*2)->second);
        else
            substitution.push_back(cudd.ReadVars(i));
    }

    return substitution;
}


vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}


void print_aiger_like_order(const Cudd& cudd) {
    string best_order = cudd.OrderString();
    auto str_vars = split(best_order, ' ');                                 MASSERT(str_vars.size() == (unsigned) cudd.ReadSize(), "");
    for (auto const & str_v : str_vars)
        cout << stoi(str_v.substr(1))*2 << " ";
    cout << endl;
}


void reorder_opt(Cudd& cudd) {                                              L_INF("reorder_opt, original size: " << cudd.ReadNodeCount());
    int count = 0;                                                          Timer timer;
    while (count != cudd.ReadNodeCount()) {
        count = (int) cudd.ReadNodeCount();
        cudd.ReduceHeap(CUDD_REORDER_SIFT_CONVERGE, 0);                     L_INF("REORDERING: SIFT_CONV: SIZE: " << cudd.ReadNodeCount());
    }                                                                       L_INF("first reordering while took (sec): " << timer.sec_restart());

    // now add annealing
    count = 0;
    while (count != cudd.ReadNodeCount()) {
        count = (int) cudd.ReadNodeCount();
        cudd.ReduceHeap(CUDD_REORDER_ANNEALING, 0);                           L_INF("REORDERING: ANNLNG: SIZE: " << cudd.ReadNodeCount());
        cudd.ReduceHeap(CUDD_REORDER_SIFT_CONVERGE, 0);                       L_INF("REORDERING: SIFT_CONV: SIZE: " << cudd.ReadNodeCount());
    }                                                                         L_INF("second reordering while took (sec): " << timer.sec_restart());
}


BDD Synth::calc_win_region() {
    /** Calculate a winning region for the safety game: win = greatest_fix_point.X [not_error & pre_sys(X)]
        :return: BDD representing the winning region
    **/

    BDD new_ = cudd.bddOne();
    for (unsigned i = 1; ; ++i) {                                             L_INF("calc_win_region: iteration " << i << "..");
        BDD curr = new_;

        new_ = pre_sys(curr);

        if (!(init <= new_))  // checking containment
            return cudd.bddZero();

        if (new_ == curr)
            return new_;

        //TODO: after the first/second iteration do expensive reordering? and build groups?
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
        aiger_symbol *aiger_input = aiger_is_input(aiger_spec, STRIP_LIT(controls[i].NodeReadIndex()*2));
        L_INF("getting output function for " << aiger_input->name);

        // the current control is [0], others = [1..]
        swap(controls[0], controls[i]);
        BDD c = controls[0];

        BDD c_arena;
        if (controls.size() >= 2) {
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

        nondet_strategy = nondet_strategy & make_bdd_eq(c, c_model);   // TODO: call substitute instead

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
unsigned Synth::walk(DdNode *a_dd) {
    // caching
    static unordered_map<DdNode*, unsigned> cache;
    {
        auto cached_lit = cache.find(Cudd_Regular(a_dd));
        if (cached_lit != cache.end())
            return Cudd_IsComplement(a_dd) ? NEGATED(cached_lit->second) : cached_lit->second;
    }
    // end of caching

    if (Cudd_IsConstant(a_dd))
        return (unsigned) (a_dd == cudd.bddOne().getNode());  // in aiger: 0 is False and 1 is True

    // get an index of variable,
    // all variables used in BDDs are also present in AIGER
    unsigned a_lit = Cudd_NodeReadIndex(a_dd)*2;

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

    unsigned and_lit = get_optimized_and_lit(n_a_t_lit, n_na_e_lit);

    unsigned res = NEGATED(and_lit);

    cache[Cudd_Regular(a_dd)] = res;  // caching

    if (Cudd_IsComplement(a_dd))
        res = NEGATED(res);

    return res;
}


/* Update aiger spec with a definition of `c_signal` */
void Synth::model_to_aiger(BDD &c_signal, BDD &func) {
    unsigned c_lit = c_signal.NodeReadIndex()*2;

    unsigned func_as_aiger_lit = walk(func.getNode());

    aiger_redefine_input_as_and(aiger_spec, c_lit, func_as_aiger_lit, func_as_aiger_lit);
}


bool first_cmp (pair<unsigned, unsigned> a, pair<unsigned, unsigned> b) { return (a.first > b.first); /* > means often-first */ }


void print_order_frequencies(Grapher& grapher, Cudd& cudd, aiger* aiger_spec) {
    string best_order = cudd.OrderString();
    auto str_vars = split(best_order, ' ');                             MASSERT(str_vars.size() == (unsigned) cudd.ReadSize(), "");
    vector<int> int_vars;
    for (auto const & str_v : str_vars)
        int_vars.push_back(stoi(str_v.substr(1)));
    cout << endl;
    cout << "frequencies of best order:" << endl;
    for (auto int_v : int_vars) {
        cout << grapher.freq_map[int_v * 2];
        auto s = aiger_is_input(aiger_spec, (unsigned int) (int_v * 2));
        if (s) {
            cout << "\t*";
            if (string(s->name).find("controllable") != string::npos)
                cout << "*";
        }
        cout << endl;
    }
    cout << endl;
}




vector<int> compute_permutation(Grapher& grapher, Cudd& cudd, aiger* aiger_spec) {
    /* Assumes the variables are already created. */

    vector<pair<unsigned, unsigned> > freq_lit_vector;
    for (unsigned i = 0; i < aiger_spec->num_inputs; ++i) {
        auto lit = aiger_spec->inputs[i].lit;
        freq_lit_vector.push_back(make_pair(grapher.freq_map[lit], lit));
    }
    for (unsigned i = 0; i < aiger_spec->num_latches; ++i) {
        auto lit = aiger_spec->latches[i].lit;
        freq_lit_vector.push_back(make_pair(grapher.freq_map[lit], lit));
    }

    sort(freq_lit_vector.begin(), freq_lit_vector.end(), first_cmp);

    MASSERT(((unsigned)cudd.ReadSize()) == aiger_spec->num_inputs + aiger_spec->num_latches+1, "should not happen");

    unordered_set<unsigned> all_var_indices;
    for (unsigned i = 0; i < (unsigned)cudd.ReadSize(); ++i)
        all_var_indices.insert(i);

    vector<int> permutation;
    for (unsigned i = 0; i < (unsigned)cudd.ReadSize(); ++i) {
        if (i < (aiger_spec->num_inputs+aiger_spec->num_latches)) {
            permutation.push_back(freq_lit_vector[i].second/2);
            all_var_indices.erase(freq_lit_vector[i].second/2);
//                cout << "perm level " << i << " lit freq " << freq_lit_vector[i].first << " setting var(lit/2) " << freq_lit_vector[i].second/2 << endl;
        }
        else {
            auto some_random_var = *all_var_indices.begin();                         L_INF("ALERT!!! who did create this var? " << some_random_var);
            all_var_indices.erase(all_var_indices.begin());
            permutation.push_back(some_random_var);
        }
    }
    MASSERT(all_var_indices.empty(), "should be empty");

    return permutation;
}


bool Synth::run() {
    L_INF("synthesize..");

    Grapher grapher(this->aiger_spec);                              timer.sec_restart();
    grapher.compute_deps();                                         L_INF("calculating deps graph took (sec): " << timer.sec_restart());
//                                               grapher.dump_dot();
//                                               cudd.EnableReorderingReporting();

    // Create all variables. _tmp ensures that BDD have positive refs.
    vector<BDD> _tmp;
    for (unsigned i=0; i < aiger_spec->num_inputs+aiger_spec->num_latches+1; ++i)
        _tmp.push_back(cudd.bddVar(i));

//    Cudd_MakeTreeNode(cudd.getManager(), 35, 5, MTR_DEFAULT);

//    vector<int> permutation = compute_permutation(grapher, cudd, aiger_spec);
//    MASSERT(permutation.size() == (unsigned) cudd.ReadSize(), "");

    /*
    L_INF("frequencies of latches");
    for (unsigned i = 0; i < aiger_spec->num_latches; ++i) {
        auto lit = aiger_spec->latches[i].lit;
        cout << "latch lit " << lit << " : " << grapher.freq_map[lit] << endl;
    }
    L_INF("frequencies of inputs");
    for (unsigned i = 0; i < aiger_spec->num_inputs; ++i) {
        auto lit = aiger_spec->inputs[i].lit;
        cout << "input lit " << lit << " : " << grapher.freq_map[lit] << endl;
    }
    */

//    cudd.ShuffleHeap(permutation.data());

    /*
    vector<BDD> nodes;
    nodes.push_back(error);
    vector<string> names;
    vector<const char*> names_;

    names.push_back(string("weird0"));
    names_.push_back(names[0].data());

    for (unsigned i = 1; i < aiger_spec->num_latches + aiger_spec->num_inputs+1; ++i) {
        names.push_back(to_string(i));
        if (aiger_is_input(aiger_spec, i*2)) {
            auto s = aiger_is_input(aiger_spec, i*2);
            if (s->name)
                names.push_back(string(s->name));
            else
                names.push_back(to_string(i*2));
        }
        if (aiger_is_latch(aiger_spec, i*2)) {
            auto s = aiger_is_latch(aiger_spec, i*2);
            if (s->name)
                names.push_back(string(s->name));
            else
                names.push_back(to_string(i*2));
        }

        names_.push_back(names[i].data());

    }
    cudd.DumpDot(nodes, names_.data(), NULL);
    cout << cudd.OrderString() << endl;
    exit(0);
    */
                                                                timer.sec_restart();
    introduce_error_bdd();                                      L_INF("introduce_error_bdd took (sec): " << timer.sec_restart());
    compose_init_state_bdd();
    compose_transition_vector();                                L_INF("calc_trans_rel took (sec): " << timer.sec_restart());

//                                               reorder_opt(cudd);
//                                               print_aiger_like_order(cudd);

                                                                timer.sec_restart();
    BDD win_region = calc_win_region();                         L_INF("calc_win_region took (sec): " << timer.sec_restart());

//                                               print_aiger_like_order(cudd);
//                                               Cudd_MakeTreeNode(cudd.getManager(), 5, 8, MTR_FIXED);

//                                               reorder_opt(cudd);
//                                               cout << "optimal order after calc_win_region" << endl;
//                                               print_aiger_like_order(cudd);
//                                                 cout << cudd.ReadNodeCount() << endl;

    if (win_region.IsZero())
        return 0;

    BDD non_det_strategy = get_nondet_strategy(win_region);

    win_region = cudd.bddZero();     // TODO: does not seem to help (also need to kill error/latches) -- transfer to new manager instead?

    vector<BDD> models = extract_output_funcs(non_det_strategy);   L_INF("extract_output_funcs took (sec): " << timer.sec_restart());
                                                                   L_INF("order_after_circuit_extraction: " << cudd.OrderString());

    vector<BDD> c_signals = get_controllable_vars_bdds();
    for (unsigned i = 0; i < models.size(); ++i)
        model_to_aiger(c_signals[i], models[i]);
                                                                   L_INF("model_to_aiger took (sec): " << timer.sec_restart());
                                                                   L_INF("circuit size: " << (aiger_spec->num_ands + aiger_spec->num_latches));
    int res = 1;
    if (output_file_name == "stdout")
        res = aiger_write_to_file(aiger_spec, aiger_ascii_mode, stdout);
    else if (!output_file_name.empty()) {                                       L_INF("writing a model to " << output_file_name);
        res = aiger_open_and_write_to_file(aiger_spec, output_file_name.c_str());
    }

    MASSERT(res, "Could not write result file");

    return 1;
}


Synth::Synth(const string &aiger_file_name,
             const string &output_file_name,
             bool calc_init_order) : output_file_name(output_file_name),
                                     calc_init_order(calc_init_order) {

    cudd.AutodynEnable(CUDD_REORDER_SIFT);

    aiger_spec = aiger_init();
    const char *err = aiger_open_and_read_from_file(aiger_spec, aiger_file_name.c_str());
    MASSERT(err == NULL, err);
}


Synth::~Synth() {
    aiger_reset(aiger_spec);
}
