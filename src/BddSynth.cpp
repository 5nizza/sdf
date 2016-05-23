#include "BddSynth.hpp"

#include "Logger.hpp"


#define IS_NEGATED(l) ((l & 1) == 1)
#define STRIP_LIT(lit) (lit & ~1)
#define NEGATED(lit) (lit ^ 1)


int BddSynth::is_fake_error_lit(unsigned lit) {
  return error_fake_latch.lit == lit;
}


BDD BddSynth::get_bdd_for_value(unsigned lit)
{
  /* - lit is variable index with sign
   */

  unsigned stripped_lit = STRIP_LIT(lit);

  // note that we faked error latch
  BDD res;

  if (stripped_lit == 0) {
    res = cudd->bddZero();
  }
  else if (is_fake_error_lit(stripped_lit) ||
           aiger_is_input(aiger_spec, stripped_lit) ||
           aiger_is_latch(aiger_spec, stripped_lit)) {
    res = cudd->bddVar(stripped_lit);  /*internal mapping of `interfaces'*/
  }
  else { // aiger_and
    aiger_and* and_ = aiger_is_and(aiger_spec, stripped_lit);
    res = get_bdd_for_value(and_->rhs0) & get_bdd_for_value(and_->rhs1);
  }

  if (IS_NEGATED(lit))
    res = !res;

  return res;
}


vector<BDD> BddSynth::get_bdd_vars(int(*filter_func)(char *)) {
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


int starts_with_controllable(char* name) { return 0 == string(name).find("controllable"); }
int not_starts_with_controllable(char* name) { return 0 != string(name).find("controllable"); }


vector<BDD> BddSynth::get_controllable_vars_bdds() {
  return get_bdd_vars(starts_with_controllable);
}


vector<BDD> BddSynth::get_uncontrollable_output_bdds() {
  return get_bdd_vars(not_starts_with_controllable);
}


void BddSynth::introduce_error_latch() {
  error_fake_latch.lit = (aiger_spec->maxvar + 1) * 2;
  error_fake_latch.name = (char *) "fake_error_latch";
  error_fake_latch.next = aiger_spec->outputs[0].lit;
}


BDD BddSynth::make_bdd_eq(BDD first, BDD second) {
  return (first & second) | (~first & ~second);
}


BDD BddSynth::compose_init_state_bdd() { // Initial state is 'all latches are zero'
  L_INF("compose_init_state_bdd..");

  BDD error_latch = get_bdd_for_value(error_fake_latch.lit);
  BDD init_state = ~error_latch;

  for(unsigned i=0; i < aiger_spec->num_latches; i++) {
    BDD latch_var = get_bdd_for_value(aiger_spec->latches[i].lit);
    init_state = init_state & ~latch_var;
  }

  return init_state;
}


BDD BddSynth::get_primed_variable_as_bdd(unsigned lit) {
  unsigned stripped_lit = STRIP_LIT(lit);
  return cudd->bddVar(stripped_lit + 1);  // in aiger odd vars are not used as names of latches/inputs
}


BDD BddSynth::compose_transition_bdd() {
  /** return: BDD representing transition function of the spec: ``T(x,i,c,x')`` **/
  L_INF("compose_transition_bdd..");

  BDD transition = cudd->bddOne();
  for (unsigned i=0; i <= aiger_spec->num_latches; ++i) {
    aiger_symbol latch = (i < aiger_spec->num_latches) ?
                         aiger_spec->latches[i] : error_fake_latch;
    BDD nextval_func = get_bdd_for_value(latch.next);
    BDD nextval_var = get_primed_variable_as_bdd(latch.lit);
    BDD latch_transition = make_bdd_eq(nextval_var, nextval_func);

    transition = transition & latch_transition;
  }

  return transition;
}


vector<BDD> BddSynth::get_all_latches_as_bdds() {
  vector<BDD> result;

  for (unsigned i = 0; i < aiger_spec->num_latches; ++i)
    result.push_back(get_bdd_for_value(aiger_spec->latches[i].lit));

  result.push_back(get_bdd_for_value(error_fake_latch.lit));

  return result;
}


BDD BddSynth::prime_unprime_latches_in_bdd(BDD bdd, int should_prime) {
  vector<BDD> cur_vars = get_all_latches_as_bdds();

  vector<BDD> primed_vars;
  for (unsigned i=0; i < cur_vars.size(); ++i)
    primed_vars.push_back(get_primed_variable_as_bdd(cur_vars[i].NodeReadIndex()));

  BDD result;
  if (should_prime)
    result = bdd.SwapVariables(cur_vars, primed_vars);
  else
    result = bdd.SwapVariables(primed_vars, cur_vars);

  return result;
}


BDD BddSynth::prime_latches_in_bdd(BDD bdd) { return prime_unprime_latches_in_bdd(bdd, 1); }
BDD BddSynth::unprime_latches_in_bdd(BDD bdd) { return prime_unprime_latches_in_bdd(bdd, 0); }


BDD BddSynth::pre_sys_bdd(BDD dst_states, BDD transition) {  //COMP: faster is possible
  /** Calculate predecessor states of given states
  :return: BDD representation of predecessor states

  :note: we need to prime dst_states before the computation (why?)
  **/

  BDD primed_dst_states = prime_latches_in_bdd(dst_states);
  BDD intersection = transition & primed_dst_states;

  BDD exist_c;
  vector<BDD> controllable = get_controllable_vars_bdds();
  if (!controllable.empty()) {
    // ∃c tau(t,i,t',c)
    BDD vars_cube = cudd->bddComputeCube(controllable.data(), NULL, (int) controllable.size());
    exist_c = intersection.ExistAbstract(vars_cube);
  }
  else {
    exist_c = intersection;
  }

  vector<BDD> all_latches = get_all_latches_as_bdds();  //ALL including fake error latch
  BDD all_latches_cube = cudd->bddComputeCube(all_latches.data(), NULL, (int)all_latches.size());
  BDD p_all_latches_cube = prime_latches_in_bdd(all_latches_cube);
  // ∃c ∃t'  tau(t,i,t',c)
  BDD exist_c_exist_next = exist_c.ExistAbstract(p_all_latches_cube);

  BDD forall_i;
  vector<BDD> uncontrollable = get_uncontrollable_output_bdds();
  if (!uncontrollable.empty()) {
    // ∀i ∃c ∃t'  tau(t,i,t',c);
    BDD uncontrollable_cube = cudd->bddComputeCube(uncontrollable.data(), NULL, (int)uncontrollable.size());
    forall_i = exist_c_exist_next.UnivAbstract(uncontrollable_cube);
  }
  else {
    forall_i = exist_c_exist_next;
  }

  return forall_i;
}


BDD BddSynth::calc_win_region(BDD init, BDD transition) {
  /** Calculate a winning region for the safety game: win = greatest_fix_point.X [not_error & pre_sys(X)]
      :return: BDD representing the winning region (win_region = win_region(x))
  **/

  L_INF("calc_win_region..");
  BDD not_err = ~get_bdd_for_value(error_fake_latch.lit);
  BDD new_ = cudd->bddOne();
  while (1) {
    BDD curr = new_;

    BDD pre_sys = pre_sys_bdd(curr, transition);
    new_ = not_err & pre_sys;

    BDD init_intersection = new_ & init;
    if (init_intersection == cudd->bddZero()) {
      return cudd->bddZero();
    }

    if (new_ == curr) {
      return new_;
    }
  }
}


BDD BddSynth::get_nondet_strategy(BDD win_region, BDD transition) {
  /** Get non-deterministic strategy from the winning region.
  If the system outputs controllable values that satisfy this non-deterministic strategy, then the system wins.
  I.e., a non-deterministic strategy describes for each state all possible plausible output values:

  ``strategy(x,i,c) = ∃x' W(x) & T(x,i,c,x') & W(x') ``

  :return: non deterministic strategy bdd
  :note: The strategy is non-deterministic -- determinization is done later.
  **/

  L_INF("get_nondet_strategy..");

  BDD primed_win_region = prime_latches_in_bdd(win_region);
  BDD intersection = primed_win_region & transition & win_region;

  vector<BDD> all_latches = get_all_latches_as_bdds();
  BDD all_latches_cube = cudd->bddComputeCube(all_latches.data(), NULL, (int) all_latches.size());
  BDD primed_all_latches_cube = prime_latches_in_bdd(all_latches_cube);
  BDD strategy = intersection.ExistAbstract(primed_all_latches_cube);

  return strategy;
}


/* The result vector respects the order of the controllable variables */
vector<BDD> BddSynth::extract_output_funcs(BDD nondet_strategy) {
  L_INF("extract_output_funcs..");

  vector<BDD> models;

  // BDD controls[aiger_spec->num_inputs];
  vector<BDD> controls = get_controllable_vars_bdds();
  for (unsigned i = 0; i < controls.size(); ++i) {
    aiger_symbol* aiger_input = aiger_is_input(aiger_spec, STRIP_LIT(controls[i].NodeReadIndex()));
    L_INF("getting output function for " << aiger_input->name << endl);

    // TODO: beautify: the current control is [0], others = [1..]
    swap(controls[0], controls[i]);
    BDD c = controls[0];

    BDD c_arena;
    if (controls.size()) {
      BDD cube = cudd->bddComputeCube(&controls[1], NULL, (int) (controls.size() - 1));
      c_arena = nondet_strategy.ExistAbstract(cube);
    }
    else { //special case of a single control
      c_arena = nondet_strategy;
    }

    BDD c_can_be_true = c_arena.Cofactor(c);
    BDD c_can_be_false = c_arena.Cofactor(~c);

    // We need to intersect with can_be_true to narrow the search.
    // Negation can cause including states from !W (with err=1)
    BDD c_must_be_true = ~c_can_be_false & c_can_be_true;
    BDD c_must_be_false = c_can_be_false & ~c_can_be_true;

    BDD c_care_set = c_must_be_true | c_must_be_false;

    // We use 'restrict' operation, but we could also do just:
    // c_model = must_be_true -> care_set
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
//    c_model.PrintMinterm();

    nondet_strategy = nondet_strategy & make_bdd_eq(c, c_model);
  }

  swap(controls[0], controls[controls.size()-1]);

  return models;
}


unsigned BddSynth::next_lit()
{
  /* return: next possible to add to the spec literal */
  return (aiger_spec->maxvar + 1) * 2;
}


unsigned BddSynth::get_optimized_and_lit(unsigned a_lit, unsigned b_lit) {
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
:warning: variables in Cudd nodes may be complemented, check with: ``node.IsComplement()``
*/
unsigned BddSynth::walk(DdNode* a_dd) {
  if (Cudd_IsConstant(a_dd))
    return (unsigned) (a_dd == cudd->bddOne().getNode());  // in aiger: 0 / 1 = False / True

  // get an index of variable,
  // all variables used in BDDs are also present in AIGER
  unsigned a_lit = Cudd_NodeReadIndex(a_dd);

  DdNode* t_bdd = Cudd_T(a_dd);
  DdNode* e_bdd = Cudd_E(a_dd);

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
void BddSynth::model_to_aiger(BDD &c_signal, BDD &func) {
  unsigned c_lit = c_signal.NodeReadIndex();

  unsigned func_as_aiger_lit = walk(func.getNode());

  aiger_redefine_input_as_and(aiger_spec, c_lit, func_as_aiger_lit, func_as_aiger_lit);
}


bool BddSynth::run() {
  L_INF("synthesize..");

  introduce_error_latch();
  BDD init_state_bdd = compose_init_state_bdd();
  BDD transition_bdd = compose_transition_bdd();


//  time_t t;
//  struct tm * now;
//
//  t = time(0);   // get time now
//  now = localtime( & t );
//  cout << now->tm_hour << ":" << now->tm_min << endl;
//
//  std::cout << "REORDERING: ORIGINAL: SIZE: " << this->cudd->ReadNodeCount() << std::endl;
//
//  int count = 0;
//  while (count != this->cudd->ReadNodeCount()) {
//    count = this->cudd->ReadNodeCount();
//
//    this->cudd->ReduceHeap(CUDD_REORDER_SIFT_CONVERGE);
//    std::cout << "REORDERING: SIFT_CONV: SIZE: " << this->cudd->ReadNodeCount() << std::endl;
//
//    this->cudd->ReduceHeap(CUDD_REORDER_LINEAR_CONVERGE);
//    std::cout << "REORDERING: LIN_CONV: SIZE: " << this->cudd->ReadNodeCount() << std::endl;
//  }
//
//  // now add annealing
//  count = 0;
//  while (count != this->cudd->ReadNodeCount()) {
//    count = this->cudd->ReadNodeCount();
//
//    // very slow, but very good results!
//    this->cudd->ReduceHeap(CUDD_REORDER_ANNEALING);
//    std::cout << "REORDERING: ANNLNG: SIZE: " << this->cudd->ReadNodeCount() << std::endl;
//
//    this->cudd->ReduceHeap(CUDD_REORDER_SIFT_CONVERGE);
//    std::cout << "REORDERING: SIFT_CONV: SIZE: " << this->cudd->ReadNodeCount() << std::endl;
//
//    this->cudd->ReduceHeap(CUDD_REORDER_LINEAR_CONVERGE);
//    std::cout << "REORDERING: LIN_CONV: SIZE: " << this->cudd->ReadNodeCount() << std::endl;
//  }
//
////  this->cudd->ReduceHeap(CUDD_REORDER_EXACT);
////  std::cout << "after6: " << this->cudd->ReadNodeCount() << std::endl;
//
//  t = time(0);   // get time now
//  now = localtime( & t );
//  cout << now->tm_hour << ":" << now->tm_min << endl;
//
//  std::cout << "VARIABLE ORDER" << std::endl;
//  std::cout << this->cudd->OrderString() << std::endl;
//  exit(0);

  BDD win_region = calc_win_region(init_state_bdd, transition_bdd);

  if (win_region.IsZero())
    return 0;

  BDD non_det_strategy = get_nondet_strategy(win_region, transition_bdd);
  vector<BDD> models = extract_output_funcs(non_det_strategy);
  vector<BDD> c_signals = get_controllable_vars_bdds();

  for (unsigned i = 0; i < models.size(); ++i)
    model_to_aiger(c_signals[i], models[i]);

  int res = (_output_file_name == "stdout") ?  //TODO: magic constant
            aiger_write_to_file(aiger_spec, aiger_ascii_mode, stdout):
            aiger_open_and_write_to_file(aiger_spec, _output_file_name.c_str());
  MASSERT(res, "Could not write result file");

  return 1;
}


BddSynth::BddSynth(const string& aiger_file_name,
                   const string& output_file_name): _output_file_name(output_file_name) {
  cudd = new Cudd();
  cudd->AutodynEnable(CUDD_REORDER_SIFT);

  aiger_spec = aiger_init();
  const char* err = aiger_open_and_read_from_file(aiger_spec, aiger_file_name.c_str());
  MASSERT(err==NULL, err);

  std::cout << "FILE NAME" << std::endl;
  std::cout << aiger_file_name << std::endl;
}


BddSynth::~BddSynth() {
  delete cudd;
  aiger_reset(aiger_spec);
}
