#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <vector>
#include <algorithm>
#include <spdlog/spdlog.h>

#include "Synth.hpp"


#define IS_NEGATED(l) ((l & 1) == 1)
#define STRIP_LIT(lit) (lit & ~1)
#define NEGATED(lit) (lit ^ 1)


#define L_INF(message) {spdlog::get("console")->info() << message;}


#define hmap unordered_map
typedef unsigned uint;
typedef vector<uint> VecUint;
typedef unordered_set<uint> SetUint;


/*
class Grapher {
public:
    Grapher() {}

    hmap<uint, SetUint> deps;  // aiger_uint_lit -> set of aiger_uint_lit // for latches -- use latch.next
    hmap<uint, uint> freq_map;

    SetUint
    _get_deps(uint lit, aiger *spec) {
        auto stripped_lit = STRIP_LIT(lit);

        if (deps.find(stripped_lit) == deps.end()) {
            auto and_ = aiger_is_and(spec, stripped_lit);                       MASSERT(and_, "impossible");
            auto rhs0_deps = _get_deps(and_->rhs0, spec);
            auto rhs1_deps = _get_deps(and_->rhs1, spec);

            SetUint a_deps;
            a_deps.insert(rhs0_deps.begin(), rhs0_deps.end());
            a_deps.insert(rhs1_deps.begin(), rhs1_deps.end());

            deps[stripped_lit] = a_deps;
        }

        return deps[stripped_lit];
    }

    void compute_deps(aiger* spec) {
        deps.clear();
        freq_map.clear();

        deps[0] = SetUint();  // value 'false' (need to insert this value, since we later call 'find')

        for (uint i=0; i < spec->num_latches; ++i)
            deps[spec->latches[i].lit].insert(spec->latches[i].lit);

        for (uint i=0; i < spec->num_inputs; ++i)
            deps[spec->inputs[i].lit].insert(spec->inputs[i].lit);

        // above we initialized values for latches and inputs ('real' variables)
        // now we calculate dependencies between them

        for (uint i=0; i < spec->num_ands; ++i)
            deps[spec->ands[i].lhs] = _get_deps(spec->ands[i].lhs, spec);

        // at this point deps of all 'real' vars are computed:
        // - input vars do not have deps,
        // - latches vars are computed because dependencies of their next literals are computed
        _compute_freq_map(spec);
    }

    void _compute_freq_map(aiger *spec) {
        for (uint i=0; i < spec->num_latches; ++i) {
            auto latch_deps = deps[STRIP_LIT(spec->latches[i].next)];
            for (auto const & elem : latch_deps)
                ++freq_map[elem];
        }

        auto output_deps = deps[STRIP_LIT(spec->outputs[0].lit)];
        for (auto const & elem : output_deps)
            ++freq_map[elem];
    }

    // TODO:
//    uint is_referenced_only_once(uint and_lit) { }
//
//     TODO:
//    void distance(uint lit1, uint lit2) { }
//
    void dump_dot(aiger* spec) {
        cout << "digraph latch_graph { rankdir=BT;" << endl;

        for (uint i=0; i < spec->num_inputs; i++)
            cout << spec->inputs[i].lit << "[shape=triangle];" << endl;
        for (uint i=0; i < spec->num_latches; i++)
            cout << spec->latches[i].lit << "[shape=diamond, color=magenta];" << endl;
        cout << STRIP_LIT(spec->outputs[0].lit) << "[shape=triangle, color=blue];" << endl;

        for (uint i=0; i < spec->num_latches; i++)
            for (auto const & src : deps[STRIP_LIT(spec->latches[i].next)])
                cout << src << "->" << spec->latches[i].lit << ";" << endl;

        for (auto const & src : deps[STRIP_LIT(spec->outputs[0].lit)])
            cout << src << "->" << STRIP_LIT(spec->outputs[0].lit) << ";" << endl;

        cout << "}" << endl;
    }
};
 */



BDD Synth::get_bdd_for_sign_lit(uint lit) {
    /* lit is an AIGER variable index with a 'sign' */

    if (bdd_by_aiger_unlit.find(STRIP_LIT(lit)) != bdd_by_aiger_unlit.end()) {
        BDD res = bdd_by_aiger_unlit[STRIP_LIT(lit)];
        if (IS_NEGATED(lit))
            return ~res;
    }

    uint stripped_lit = STRIP_LIT(lit);
    BDD res;

    if (stripped_lit == 0) {
        res = cudd.bddZero();
    }
    else if (aiger_is_input(aiger_spec, stripped_lit) ||
             aiger_is_latch(aiger_spec, stripped_lit)) {
        res = cudd.ReadVars(cudd_by_aiger[stripped_lit]);
//        MASSERT(res.NodeReadIndex() == stripped_lit/2, "that bug again: impossible: " << res.NodeReadIndex() << " vs " << stripped_lit/2 );
    }
    else { // aiger_and
        aiger_and *and_ = aiger_is_and(aiger_spec, stripped_lit);
        res = get_bdd_for_sign_lit(and_->rhs0) & get_bdd_for_sign_lit(and_->rhs1);
    }

    bdd_by_aiger_unlit[stripped_lit] = res;

    return IS_NEGATED(lit) ? (~res):res;
}


vector<BDD> Synth::get_bdd_vars(bool(*filter_func)(char *)) {
    vector<BDD> result;
    for (uint i = 0; i < aiger_spec->num_inputs; ++i) {
        aiger_symbol symbol = aiger_spec->inputs[i];
        if ((*filter_func)(symbol.name)) {
            BDD out_var_bdd = get_bdd_for_sign_lit(symbol.lit);
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


vector<BDD> Synth::get_uncontrollable_vars_bdds() {
    return get_bdd_vars(not_starts_with_controllable);
}


void Synth::introduce_error_bdd() {
    error = get_bdd_for_sign_lit(aiger_spec->outputs[0].lit);
}


void Synth::compose_init_state_bdd() { // Initial state is 'all latches are zero'
    L_INF("compose_init_state_bdd..");

    init = cudd.bddOne();

    for (uint i = 0; i < aiger_spec->num_latches; i++) {
        BDD latch_var = get_bdd_for_sign_lit(aiger_spec->latches[i].lit);
        init = init & ~latch_var;
    }
}


void Synth::compose_transition_vector() {
    L_INF("compose_transition_vector..");

    for (uint i = 0; i < aiger_spec->num_latches; ++i)
        transition_func[aiger_spec->latches[i].lit] = get_bdd_for_sign_lit(aiger_spec->latches[i].next);
}


vector<BDD> Synth::get_substitution() {
    vector<BDD> substitution;

    for (uint i = 0; i < (uint)cudd.ReadSize(); ++i) {
        if (aiger_is_latch(aiger_spec, aiger_by_cudd[i]))
            substitution.push_back(transition_func.find(aiger_by_cudd[i])->second);
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

/*
void print_aiger_like_order(const Cudd& cudd) {
    string best_order = cudd.OrderString();
    auto str_vars = split(best_order, ' ');                                 MASSERT(str_vars.size() == (uint) cudd.ReadSize(), "");
    for (auto const & str_v : str_vars)
        cout << stoi(str_v.substr(1))*2 << " ";
    cout << endl;
}
 */

/*
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

*/
VecUint get_order(Cudd &cudd) {
    string order_str = cudd.OrderString();
    auto str_vars = split(order_str, ' ');                                MASSERT(str_vars.size() == (uint) cudd.ReadSize(), "");
    VecUint order;
    for (auto const & str_v : str_vars) {
        order.push_back((uint) stoi(str_v.substr(1)));
    }
    return order;
}

//string string_vector(const VecUint& v) {
//    stringstream ss;
//    for(size_t i = 0; i < v.size(); ++i)
//    {
//        if(i != 0)
//            ss << ",";
//        ss << v[i];
//    }
//    return ss.str();
//}
//
string string_set(const SetUint& v) {
    stringstream ss;
    for(auto const el : v)
        ss << el << ",";
    return ss.str();
}

template<typename Container>
struct container_hash {
    std::size_t operator()(Container const & v) const {
        size_t hash = 0;
        for (auto const el: v)
            hash ^= el;
        return hash;
    }
};


bool do_intersect(const SetUint& s1, const SetUint& s2) {
    for (auto el1 : s1)
        if (s2.find(el1) != s2.end())
            return true;
    return false;
}

//SetUint union_sets(const SetUint& s1, const SetUint& s2) {
//    SetUint merged(s1.begin(), s1.end());
//    merged.insert(s2.begin(), s2.end());
//    return merged;
//}


vector<SetUint>
get_group_candidates(const vector<VecUint >& orders,
                     uint window_size)
{
    hmap<SetUint, uint, container_hash<SetUint>>
            group_freq;

    for (auto const & order : orders) {
        for (uint idx=0; idx < order.size() - window_size; ++idx) {
            SetUint sub_group(order.begin() + idx,
                              order.begin() + idx + window_size);
            ++group_freq[sub_group];
        }
    }

    vector<SetUint> candidates;
    for (auto const & it: group_freq) {
        if (((float)it.second/orders.size()) >= 0.8)  // appears 'often'
            candidates.push_back(it.first);
    }
    return candidates;
}


void remove_intersecting(SetUint group,
                         vector<SetUint> &groups)
{
    vector<SetUint >::iterator it = groups.begin();
    while (it != groups.end()) {
        if (do_intersect(group, *it))
            it = groups.erase(it);
        else
            ++it;
    }
}


uint get_var_of_min_order_position(Cudd& cudd, const SetUint& group)
{
    uint min_var = *group.begin();
    uint min_pos = (uint)cudd.ReadPerm(min_var);
    for (auto const var : group)
        if ((uint )cudd.ReadPerm(var) < min_pos) {
            min_var = var;
            min_pos = (uint)cudd.ReadPerm(var);
        }
    return min_var;
}


void introduce_group_into_cudd(Cudd &cudd, const SetUint& group)
{
    L_INF("adding variable group to cudd: " << string_set(group));
    auto first_var_pos = get_var_of_min_order_position(cudd, group);
    cudd.MakeTreeNode(first_var_pos, (uint) group.size(), MTR_FIXED);
}


void _do_grouping(Cudd &cudd,
                  hmap<uint, vector<SetUint>> &groups_by_length,  // we modify its values
                  uint cur_group_length,
                  const VecUint& cur_order)
{
    L_INF("fixing groups of size " << cur_group_length << ". The number of groups = " << groups_by_length[cur_group_length].size());

    auto cur_groups = groups_by_length[cur_group_length];

    for (uint i = 0; i+cur_group_length < cur_order.size(); ++i) {
        SetUint candidate;
        for (uint j = 0; j < cur_group_length; ++j)
            candidate.insert(cur_order[i+j]);

        if (find(cur_groups.begin(), cur_groups.end(), candidate) != cur_groups.end()) {
            for (uint l = 2; l < cur_group_length; ++l) {
//                cout << "rm intersections " << string_set(candidate) << endl;
//                cout << "before: " << endl;
//                for (auto const v : groups_by_length[l])
//                    cout << string_set(v) << endl;
                remove_intersecting(candidate, groups_by_length[l]);  //remove from smaller groups
//                cout << "after: " << endl;
//                for (auto const v : groups_by_length[l])
//                    cout << string_set(v) << endl;
            }
            introduce_group_into_cudd(cudd, candidate);
        }
    }
}


void do_grouping(Cudd& cudd,
                 const vector<VecUint>& orders)
{
    L_INF("trying to group vars..");

    if (orders[0].size() < 5)  // window size is too big
        return;

    hmap<uint, vector<SetUint>> groups_by_length;
    groups_by_length[2] = get_group_candidates(orders, 2);
    groups_by_length[3] = get_group_candidates(orders, 3);
    groups_by_length[4] = get_group_candidates(orders, 4);
    groups_by_length[5] = get_group_candidates(orders, 5);

    L_INF("# of group candidates: of size 2 -- " << groups_by_length[2].size());
    for (auto const& g : groups_by_length[2]) {
        L_INF(string_set(g));
    }
    L_INF("# of group candidates: of size 3 -- " << groups_by_length[3].size());
    for (auto const& g : groups_by_length[3]) {
        L_INF(string_set(g));
    }
    L_INF("# of group candidates: of size 4 -- " << groups_by_length[4].size());
    L_INF("# of group candidates: of size 5 -- " << groups_by_length[5].size());

    auto cur_order = orders.back();    // we fix only groups present in the current order (because that is easier to implement)

    for (uint i = 5; i>=2; --i)  // decreasing order!
        if (!groups_by_length[i].empty())
            _do_grouping(cudd, groups_by_length, i, cur_order);
}


void update_order_if(Cudd& cudd, vector<VecUint >& orders)
{
    static uint last_nof_orderings = 0;

    if (last_nof_orderings != cudd.ReadReorderings())
        orders.push_back(get_order(cudd));

    last_nof_orderings = cudd.ReadReorderings();
}


BDD Synth::pre_sys(BDD dst) {
    /**
    Calculate predecessor states of given states.

        ∀u ∃c ∃t': tau(t,u,c,t') & dst(t') & ~error(t,u,c)

    We use the direct substitution optimization (since t' <-> BDD(t,u,c)), thus:

        ∀u ∃c: (!error(t,u,c)  &  (dst(t)[t <- bdd_next_t(t,u,c)]))

    or for Moore machines:

        ∃c ∀u: (!error(t,u,c)  &  (dst(t)[t <- bdd_next_t(t,u,c)]))

    Note that we do not replace t variables in the error bdd.

    :return: BDD of the predecessor states
    **/

    // NOTE: I tried considering two special cases: error(t,u,c) and error(t),
    //       and move error(t) outside of quantification ∀u ∃c.
    //       It slowed down..
    // TODO: try again: on the driver example

    static vector<VecUint> orders;
    static bool did_grouping = false;

    if (!did_grouping && timer.sec_from_origin() > time_limit_sec/4) { // at 0.25*time_limit we fix the order
        do_grouping(cudd, orders);
        did_grouping = true;
    }

    dst = dst.VectorCompose(get_substitution());                                update_order_if(cudd, orders);

    if (is_moore) {
        BDD result = dst.And(~error);

        vector<BDD> uncontrollable = get_uncontrollable_vars_bdds();
        if (!uncontrollable.empty()) {
            BDD uncontrollable_cube = cudd.bddComputeCube(uncontrollable.data(),
                                                          NULL,
                                                          (int)uncontrollable.size());
            result = result.UnivAbstract(uncontrollable_cube);                  update_order_if(cudd, orders);
        }

        // ∃c ∀u  (...)
        vector<BDD> controllable = get_controllable_vars_bdds();
        BDD controllable_cube = cudd.bddComputeCube(controllable.data(),
                                                    NULL,
                                                    (int)controllable.size());
        result = result.ExistAbstract(controllable_cube);                       update_order_if(cudd, orders);
        return result;
    }

    // the case of Mealy machines
    vector<BDD> controllable = get_controllable_vars_bdds();
    BDD controllable_cube = cudd.bddComputeCube(controllable.data(),
                                                NULL,
                                                (int)controllable.size());
    BDD result = dst.AndAbstract(~error, controllable_cube);                   update_order_if(cudd, orders);

    vector<BDD> uncontrollable = get_uncontrollable_vars_bdds();
    if (!uncontrollable.empty()) {
        // ∀u ∃c (...)
        BDD uncontrollable_cube = cudd.bddComputeCube(uncontrollable.data(),
                                                      NULL,
                                                      (int)uncontrollable.size());
        result = result.UnivAbstract(uncontrollable_cube);                    update_order_if(cudd, orders);
    }
    return result;
}


//int post_order_hook(DdManager* cudd, const char*, void*)
//{
//    L_INF("cudd initiated reordering: ");
//    return 1;
//}


BDD Synth::calc_win_region() {
    /** Calculate a winning region for the safety game: win = greatest_fix_point.X [not_error & pre_sys(X)]
        :return: BDD representing the winning region
    **/

    BDD new_ = cudd.bddOne();
    for (uint i = 1; ; ++i) {                                             L_INF("calc_win_region: iteration " << i << ": node count " << cudd.ReadNodeCount());
        BDD curr = new_;

        new_ = pre_sys(curr);

        if ((init & new_) == cudd.bddZero())
            return cudd.bddZero();

        if (new_ == curr)
            return new_;
    }
}


BDD Synth::get_nondet_strategy() {
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

    // TODO: do we need win_region?
    return ~error & win_region & win_region.VectorCompose(get_substitution());
}


hmap<uint,BDD> Synth::extract_output_funcs() {
    /** The result vector respects the order of the controllable variables **/

    L_INF("extract_output_funcs..");

    cudd.FreeTree();    // ordering that worked for win region computation might not work here

    hmap<uint,BDD> model_by_cuddidx;

    vector<BDD> controls = get_controllable_vars_bdds();

    while (!controls.empty()) {
        BDD c = controls.back(); controls.pop_back();

        aiger_symbol *aiger_input = aiger_is_input(aiger_spec, aiger_by_cudd[c.NodeReadIndex()]);
        L_INF("getting output function for " << aiger_input->name);

        BDD c_arena;
        if (controls.size() > 0) {
            BDD cube = cudd.bddComputeCube(controls.data(), NULL, (int)controls.size());
            c_arena = non_det_strategy.ExistAbstract(cube);
        }
        else { //no other signals left
            c_arena = non_det_strategy;
        }
        // Now we have: c_arena(t,u,c) = ∃c_others: nondet(t,u,c)
        // (i.e., c_arena talks about this particular c, about t and u)

        BDD c_can_be_true = c_arena.Cofactor(c);
        BDD c_can_be_false = c_arena.Cofactor(~c);

        BDD c_must_be_true = ~c_can_be_false & c_can_be_true;
        BDD c_must_be_false = c_can_be_false & ~c_can_be_true;
        // Note that we cannot use `c_must_be_true = ~c_can_be_false`,
        // since the negation can cause including tuples (t,i,o) that violate non_det_strategy.

        auto support_indices = cudd.SupportIndices(vector<BDD>({c_must_be_false, c_must_be_true}));
        for (auto const var_cudd_idx : support_indices) {
            auto v = cudd.ReadVars(var_cudd_idx);
            auto new_c_must_be_false = c_must_be_false.ExistAbstract(v);
            auto new_c_must_be_true = c_must_be_true.ExistAbstract(v);

            if ((new_c_must_be_false & new_c_must_be_true) == cudd.bddZero()) {
                c_must_be_false = new_c_must_be_false;
                c_must_be_true = new_c_must_be_true;
            }
        }

        // We use 'restrict' operation, but we could also just do:
        //     c_model = care_set -> must_be_true
        // but this is (presumably) less efficient (in time? in size?).
        // (intuitively, because we always set c_model to 1 if !care_set, but we could set it to 0)
        //
        // The result of restrict operation satisfies:
        //     on c_care_set: c_must_be_true <-> must_be_true.Restrict(c_care_set)

        BDD c_model = c_must_be_true.Restrict(c_must_be_true | c_must_be_false);

        model_by_cuddidx[c.NodeReadIndex()] = c_model;

        //killing node refs
        c_must_be_false = c_must_be_true = c_can_be_false = c_can_be_true = c_arena = cudd.bddZero();

        //TODO: ak: strange -- the python version for the example amba_02_9n produces a smaller circuit (~5-10 times)!
        non_det_strategy = non_det_strategy.Compose(c_model, c.NodeReadIndex());
        //non_det_strategy = non_det_strategy & ((c & c_model) | (~c & ~c_model));
    }

    return model_by_cuddidx;
}


uint Synth::next_lit() {
    /* return: next possible to add to the spec literal */
    return (aiger_spec->maxvar + 1) * 2;
}


uint Synth::get_optimized_and_lit(uint a_lit, uint b_lit) {
    if (a_lit == 0 || b_lit == 0)
        return 0;

    if (a_lit == 1 && b_lit == 1)
        return 1;

    if (a_lit == 1)
        return b_lit;

    if (b_lit == 1)
        return a_lit;

    if (a_lit > 1 && b_lit > 1) {
        uint a_b_lit = next_lit();
        aiger_add_and(aiger_spec, a_b_lit, a_lit, b_lit);
        return a_b_lit;
    }

    MASSERT(0, "impossible");
}


uint Synth::walk(DdNode *a_dd) {
    /**
    Walk given DdNode node (recursively).
    If a given node requires intermediate AND gates for its representation, the function adds them.
        Literal representing given input node is `not` added to the spec.

    :returns: literal representing input node
    **/

    // caching
    static hmap<DdNode*, uint> cache;
    {
        auto cached_lit = cache.find(Cudd_Regular(a_dd));
        if (cached_lit != cache.end())
            return Cudd_IsComplement(a_dd) ? NEGATED(cached_lit->second) : cached_lit->second;
    }
    // end of caching

    if (Cudd_IsConstant(a_dd))
        return (uint) (a_dd == cudd.bddOne().getNode());  // in aiger: 0 is False and 1 is True

    // get an index of the variable
    uint a_lit = aiger_by_cudd[Cudd_NodeReadIndex(a_dd)];

    DdNode *t_bdd = Cudd_T(a_dd);
    DdNode *e_bdd = Cudd_E(a_dd);

    uint t_lit = walk(t_bdd);
    uint e_lit = walk(e_bdd);

    // ite(a_bdd, then_bdd, else_bdd)
    // = a*then + !a*else
    // = !(!(a*then) * !(!a*else))
    // -> in general case we need 3 more ANDs

    uint a_t_lit = get_optimized_and_lit(a_lit, t_lit);

    uint na_e_lit = get_optimized_and_lit(NEGATED(a_lit), e_lit);

    uint n_a_t_lit = NEGATED(a_t_lit);
    uint n_na_e_lit = NEGATED(na_e_lit);

    uint and_lit = get_optimized_and_lit(n_a_t_lit, n_na_e_lit);

    uint res = NEGATED(and_lit);

    cache[Cudd_Regular(a_dd)] = res;

    if (Cudd_IsComplement(a_dd))
        res = NEGATED(res);

    return res;
}


void Synth::model_to_aiger(const BDD &c_signal, const BDD &func) {
    /// Update AIGER spec with a definition of `c_signal`

    uint c_lit = aiger_by_cudd[c_signal.NodeReadIndex()];
    string output_name = string(aiger_is_input(aiger_spec, c_lit)->name);  // save the name before it is freed

    uint func_as_aiger_lit = walk(func.getNode());

    aiger_redefine_input_as_and(aiger_spec, c_lit, func_as_aiger_lit, func_as_aiger_lit);

    if (print_full_model)
        aiger_add_output(aiger_spec, c_lit, output_name.c_str());
}


//bool first_cmp (pair<uint, uint> a, pair<uint, uint> b) { return (a.first > b.first); /* > means often-first */ }

/*
void print_order_frequencies(Grapher& grapher, Cudd& cudd, aiger* aiger_spec) {
    string best_order = cudd.OrderString();
    auto str_vars = split(best_order, ' ');                             MASSERT(str_vars.size() == (uint) cudd.ReadSize(), "");
    vector<int> int_vars;
    for (auto const & str_v : str_vars)
        int_vars.push_back(stoi(str_v.substr(1)));
    cout << endl;
    cout << "frequencies of best order:" << endl;
    for (auto int_v : int_vars) {
        cout << grapher.freq_map[int_v * 2];
        auto s = aiger_is_input(aiger_spec, (uint) (int_v * 2));
        if (s) {
            cout << "\t*";
            if (string(s->name).find("controllable") != string::npos)
                cout << "*";
        }
        cout << endl;
    }
    cout << endl;
}
*/

/*
vector<int> compute_permutation(Grapher& grapher, Cudd& cudd, aiger* aiger_spec) {
    // Assumes the variables are already created.

    vector<pair<uint, uint>> freq_lit_vector;
    for (uint i = 0; i < aiger_spec->num_inputs; ++i) {
        auto lit = aiger_spec->inputs[i].lit;
        freq_lit_vector.push_back(make_pair(grapher.freq_map[lit], lit));
    }
    for (uint i = 0; i < aiger_spec->num_latches; ++i) {
        auto lit = aiger_spec->latches[i].lit;
        freq_lit_vector.push_back(make_pair(grapher.freq_map[lit], lit));
    }

    sort(freq_lit_vector.begin(), freq_lit_vector.end(), first_cmp);

    MASSERT(((uint)cudd.ReadSize()) == aiger_spec->num_inputs + aiger_spec->num_latches+1, "should not happen");

    SetUint all_var_indices;
    for (uint i = 0; i < (uint)cudd.ReadSize(); ++i)
        all_var_indices.insert(i);

    vector<int> permutation;
    for (uint i = 0; i < (uint)cudd.ReadSize(); ++i) {
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


void print_set(SetUint& s, aiger* spec) {
    for (auto const e: s) {
        cout << e;
        auto latch = aiger_is_latch(spec, e);
        if (latch) {
            cout << "(" << latch->next << ")";
        }
        cout << endl;

    }
    cout << endl;
}
*/

class Cleaner {
public:
    aiger* spec;
    Cleaner(aiger* spec): spec(spec) { }
    ~Cleaner() { aiger_reset(spec); }
};


void init_cudd(Cudd& cudd)
{
//    cudd.Srandom(827464282);  // for reproducibility (?)
    cudd.AutodynEnable(CUDD_REORDER_SIFT);
//    cudd.EnableReorderingReporting();
}


/* something does not work here
BDD compute_reachable(aiger* spec, BDD init_orig, hmap<unsigned,BDD>& transition_func_orig, BDD error_orig,
                      Cudd& cudd_orig)
{

#define AIG_2_PRIM_IDX(aiger_lit)  (aiger_lit/2 + spec->num_latches)
#define AIG_2_CUR_IDX(aiger_lit)  (aiger_lit/2)
#define PRIM_IDX_2_AIG_LIT(aiger_lit)  (aiger_lit/2 + spec->num_inputs + spec->num_latches)
#define CUR_IDX_2_AIG(aiger_lit)  (aiger_lit/2)

    L_INF("compute_reachable...");
    Cudd cudd2;
    init_cudd(cudd2);

    vector<BDD> _tmp;
    for (uint i=0; i < spec->num_inputs + 2*spec->num_latches + 1; ++i)
        _tmp.push_back(cudd2.bddVar(i));

    VecUint permutation;
    permutation.push_back(0);

    for (uint i=0; i < spec->num_latches; ++i) {
        permutation.push_back(AIG_2_CUR_IDX(spec->latches[i].lit));
        permutation.push_back(AIG_2_PRIM_IDX(spec->latches[i].lit));
    }

    for (uint i=0; i < spec->num_inputs; ++i)
        permutation.push_back(AIG_2_CUR_IDX(spec->inputs[i].lit));

    for (auto c : permutation)
        cout << c << endl;
    cudd2.ShuffleHeap((int*)permutation.data());

    for (uint i=1; i < spec->num_latches + 1; ++i)
        cudd2.MakeTreeNode(i*2-1, 2, MTR_FIXED);

    auto init = init_orig.Transfer(cudd2);

    cout << "transferring error.." << endl;
    auto error = error_orig.Transfer(cudd2);

    cout << "computing trans rel .. " << endl;
    auto trans_rel = cudd2.bddOne();
    for (auto const& it : transition_func_orig) {
        cout << "adding & for var" << it.first << endl;
        auto newBdd = it.second.Transfer(cudd2);
        trans_rel &= ~(cudd2.ReadVars(AIG_2_PRIM_IDX(it.first)).Xor(newBdd));  //aka: l' <-> bdd
    }

    vector<BDD> all_inputs;
    for (uint i = 0; i < spec->num_inputs; ++i) {
        auto aig_inp = spec->inputs[i];
        all_inputs.push_back(cudd2.ReadVars(AIG_2_CUR_IDX(aig_inp.lit)));
    }
    auto all_inputs_cube = cudd2.bddComputeCube(all_inputs.data(), NULL, (int) all_inputs.size());

    vector<BDD> all_cur_latches;
    vector<BDD> all_prim_latches;
    for (uint i = 0; i < spec->num_latches; ++i) {
        auto aig_l = spec->latches[i];
        all_cur_latches.push_back(cudd2.ReadVars(AIG_2_CUR_IDX(aig_l.lit)));
        all_prim_latches.push_back(cudd2.ReadVars(AIG_2_PRIM_IDX(aig_l.lit)));  //i.e., normal idx multiplied by 2
    }
    auto all_cur_latches_cube = cudd2.bddComputeCube(all_cur_latches.data(), NULL, (int) all_cur_latches.size());

    auto trans_not_err = trans_rel & ~error;
    BDD curr;
    BDD new_ = init;
    do {
        curr = new_;
        cout << "another iteration.. is curr 1?" << (new_==cudd2.bddOne()) << endl;
        // new_(l') = \exists u,c,l: curr(l) & trans(l,u,c,l') & ~error(l,u,c)
        // new_(l') = \exists l: (curr(l) & (\exists u,c: (trans(l,u,c,l') & ~error(l,u,c))) )
//        new_ = curr.AndAbstract(trans_rel.AndAbstract(~error, all_inputs_cube),
//                                all_cur_latches_cube)
//               .SwapVariables(all_prim_latches, all_cur_latches);
        new_ = trans_not_err.AndAbstract(curr, all_cur_latches_cube & all_inputs_cube);
        new_ = new_.SwapVariables(all_prim_latches, all_cur_latches);

//        new_ = curr.AndAbstract(trans_rel.ExistAbstract(all_inputs_cube),
//                                all_cur_latches_cube)
//               .SwapVariables(all_prim_latches, all_cur_latches);
    } while (new_ != curr);

    L_INF("finished computing reachable states");
    cout << (new_==cudd2.bddOne()) << ", size: " << new_.nodeCount() << endl;
    L_INF("transferring it back..");
    auto res = new_.Transfer(cudd_orig);
    L_INF("finished transferring back!");
    return res;
}
 */


bool Synth::run() {
    init_cudd(cudd);

    aiger_spec = aiger_init();
    const char *err = aiger_open_and_read_from_file(aiger_spec, aiger_file_name.c_str());
    MASSERT(err == NULL, err);
    Cleaner cleaner(aiger_spec);

    // main part
    L_INF("synthesize.. number of vars = " << aiger_spec->num_inputs + aiger_spec->num_latches);
//    grapher = new Grapher();                                        timer.sec_restart();
//    grapher->compute_deps(aiger_spec);                              L_INF("calculating deps graph took (sec): " << timer.sec_restart());

    //grapher.dump_dot();
    //print_set(grapher.deps[STRIP_LIT(aiger_spec->outputs[0].lit)], aiger_spec);

    // Create all variables. _tmp ensures that BDD have positive refs.
    vector<BDD> _tmp;
    for (uint i = 0; i < aiger_spec->num_inputs + aiger_spec->num_latches; ++i)
        _tmp.push_back(cudd.bddVar(i));

    for (uint i = 0; i < aiger_spec->num_inputs; ++i) {
        auto aiger_strip_lit = aiger_spec->inputs[i].lit;
        cudd_by_aiger[aiger_strip_lit] = i;
        aiger_by_cudd[i] = aiger_strip_lit;
    }
    for (uint i = 0; i < aiger_spec->num_latches; ++i) {
        auto aiger_strip_lit = aiger_spec->latches[i].lit;
        auto cudd_idx = i + aiger_spec->num_inputs;
        cudd_by_aiger[aiger_strip_lit] = cudd_idx;
        aiger_by_cudd[cudd_idx] = aiger_strip_lit;
    }

//    vector<int> permutation = compute_permutation(grapher, cudd, aiger_spec);
//    MASSERT(permutation.size() == (uint) cudd.ReadSize(), "");

    /*
    L_INF("frequencies of latches");
    for (uint i = 0; i < aiger_spec->num_latches; ++i) {
        auto lit = aiger_spec->latches[i].lit;
        cout << "latch lit " << lit << " : " << grapher.freq_map[lit] << endl;
    }
    L_INF("frequencies of inputs");
    for (uint i = 0; i < aiger_spec->num_inputs; ++i) {
        auto lit = aiger_spec->inputs[i].lit;
        cout << "input lit " << lit << " : " << grapher.freq_map[lit] << endl;
    }
    */

    /*
    vector<BDD> nodes;
    nodes.push_back(error);
    vector<string> names;
    vector<const char*> names_;

    names.push_back(string("weird0"));
    names_.push_back(names[0].data());

    for (uint i = 1; i < aiger_spec->num_latches + aiger_spec->num_inputs+1; ++i) {
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

    compose_init_state_bdd();
    timer.sec_restart();
    compose_transition_vector();
    L_INF("calc_trans_rel took (sec): " << timer.sec_restart());
    introduce_error_bdd();
    L_INF("introduce_error_bdd took (sec): " << timer.sec_restart());

//    cout << "before comput: nof_vars = " << cudd.ReadSize() << endl;
//    reachable = compute_reachable(aiger_spec, init, transition_func, error, cudd);
//    cout << "after comput: nof_vars = " << cudd.ReadSize() << endl;


    // no need for cache
    bdd_by_aiger_unlit.clear();

//                                               reorder_opt(cudd);
//                                               print_aiger_like_order(cudd);

    timer.sec_restart();
    win_region = calc_win_region();
    L_INF("calc_win_region took (sec): " << timer.sec_restart());

//                                               print_aiger_like_order(cudd);
//                                               Cudd_MakeTreeNode(cudd.getManager(), 5, 8, MTR_FIXED);

//                                               reorder_opt(cudd);
//                                               cout << "optimal order after calc_win_region" << endl;
//                                               print_aiger_like_order(cudd);
//                                                 cout << cudd.ReadNodeCount() << endl;

    if (win_region.IsZero()) {
        cout << "UNREALIZABLE" << endl;
        return 0;
    }

    cout << "REALIZABLE" << endl;

    non_det_strategy = get_nondet_strategy();

    //cleaning non-used bdds
    win_region = cudd.bddZero();
    transition_func.clear();
    init = cudd.bddZero();
    error = cudd.bddZero();
    //

    // TODO: set time limit on reordering? or even disable it if no time?
    hmap<uint, BDD> model_by_cuddidx = extract_output_funcs();                   L_INF("extract_output_funcs took (sec): " << timer.sec_restart());

    //cleaning non-used bdds
    non_det_strategy = cudd.bddZero();
    //

    auto elapsed_sec = time_limit_sec - timer.sec_from_origin();
    if (elapsed_sec > 100) {    // leave 100sec just in case
        auto spare_time_sec = elapsed_sec - 100;
        cudd.ResetStartTime();
        cudd.IncreaseTimeLimit((unsigned long) (spare_time_sec * 1000));
        cudd.ReduceHeap(CUDD_REORDER_SIFT_CONVERGE);
        cudd.UnsetTimeLimit();
        cudd.AutodynDisable();  // just in case -- cudd hangs on timeout
    }

    for (auto const it : model_by_cuddidx)
        model_to_aiger(cudd.ReadVars((int)it.first), it.second);
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

Synth::~Synth() {
//    if (grapher)
//        delete(grapher);
//    grapher = NULL;
}
