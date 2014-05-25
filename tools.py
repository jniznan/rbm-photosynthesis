from pysb import Model, Rule, bng
from collections import defaultdict, Counter
from operator import add


def search_model(model, sort_key=None, up_to_cee=False):
    '''
    Yields all descsendants of a model. They are created by applying syntactic
    operations.

    Parameters
    ----------
    model: pysb.Model
    sort_key: function
        Function to be used for sorting the rules.
    up_to_cee : bool
        If True, then it performes only bidirectional merging and context
        enumeration elimination.
    '''
    for m in merge_bidirectional(model, sort_key):
        yield m
    for m in context_enumeration_elimination(model, sort_key):
        yield m
    if not up_to_cee:
        for m in context_elimination(model, sort_key):
            yield m
        for m in rules_removal(model):
            yield m


def merge_bidirectional(model, sort_key=None):
    '''
    Yields models that are in relation bidirectional grouping
    with the model passed in.

    Parameters
    ----------
    model: pysb.Model
    sort_key: function
        Function to be used for sorting the rules.
    '''
    lh_map = defaultdict(list)
    rh_map = defaultdict(list)
    for rule in model.rules:
        if not rule.rule_expression.is_reversible:
            lh, rh = _get_lh_rh(rule)
            lh_map[(repr(lh), repr(rh))].append(rule)
            rh_map[(repr(rh), repr(lh))].append(rule)
    for rep in set(lh_map.keys()) & set(rh_map.keys()):
        for fwd in lh_map[rep]:
            for bwd in rh_map[rep]:
                m = copy_no_rules(model)
                rules = [r for r in model.rules if r not in [fwd, bwd]]
                lh, rh = _get_lh_rh(fwd)
                lh, rh = reduce(add, lh), reduce(add, rh)
                new_rule = Rule(new_rule_name(fwd, bwd), lh <> rh,
                                fwd.rate_forward, bwd.rate_forward,
                                _export=False)
                rules.append(new_rule)
                rules.sort(key=sort_key)
                for r in rules:
                    m.add_component(r)
                yield 'bi: %s' % new_rule.name, m


def new_rule_name(*rules):
    '''
    Returns a new name for a set of rules to be merged together.
    Tries to find the largest prefix and sets it as a prefix of the new name
    in order to avoid repetition in the name.
    '''
    names = [rule.name for rule in rules]
    n = min(map(len, names))
    prefix = ''
    for i in range(n):
        chars = map(lambda n: n[i], names)
        if any(chars[0] != ch for ch in chars):
            break
        prefix += chars[0]
    ret = '__'.join(sorted(map(lambda n: n[i:].strip('_'), names)))
    if prefix:
        ret = prefix.strip('_') + '___' + ret
    return ret


def context_enumeration_elimination(model, sort_key=None):
    '''
    Yields models that are in relation context enumeration elimination
    with the model passed in.

    Parameters
    ----------
    model: pysb.Model
    sort_key: function
        Function to be used for sorting the rules.
    '''
    grouped = _group_rules_by_context_and_changes(model.rules)
    for key, group in grouped.items():
        sites = [(i, s) for i, (lname, lsites, rname, rdiffs) in enumerate(key)
                 if lname == rname for s in lsites if s not in zip(*rdiffs)[0]]
        for i, site in sites:
            # local map states of the `site` into sets of occuring contexts
            local = defaultdict(set)
            # conmap maps contexts to lists of rules that contain that context
            conmap = defaultdict(list)
            possible_states = None
            for rule in group:
                lh, _ = _get_lh_rh(rule)
                context = [(j, k, v) for j, l in enumerate(lh)
                           for k, v in l.site_conditions.items()
                           if (j, k) in sites and (k != site or i != j)]
                context = tuple(sorted(context))
                local[lh[i].site_conditions[site]].add(context)
                conmap[context].append(rule)
                if not possible_states:
                    possible_states = lh[i].monomer.site_states[site]
            ref = local[local.keys()[0]]
            if all(local[st] == ref for st in possible_states):
                #  it is enumerating, we can reduce :)
                m = copy_no_rules(model)
                rules = [r for r in model.rules if r not in group]
                new_rules = []
                for context in ref:
                    name = new_rule_name(*conmap[context])
                    rule = conmap[context][0]
                    rexp = rule.rule_expression
                    # now set the site state to None
                    lh = map(lambda p: p.copy(),
                             rexp.reactant_pattern.complex_patterns)
                    rh = map(lambda p: p.copy(),
                             rexp.product_pattern.complex_patterns)
                    lh[i].monomer_patterns[0].site_conditions[site] = None
                    rh[i].monomer_patterns[0].site_conditions[site] = None
                    rule_new = inherit_rule(rule, reduce(add, lh),
                                            reduce(add, rh), name=name)
                    new_rules.append(rule_new)
                rules.extend(new_rules)
                rules.sort(key=sort_key)
                for r in rules:
                    m.add_component(r)
                yield 'CEE: %s %s; %s' % (
                    i, site,
                    ' '.join(map(lambda r: r.name, new_rules))), m


def _group_rules_by_context_and_changes(rules):
    # assumes reactants and products are sorted the same way
    grouped = defaultdict(list)
    for rule in rules:
        lh, rh = _get_lh_rh(rule)
        n = min(len(lh), len(rh))
        key = []
        for l, r in zip(lh[:n], rh[:n]):
            lmon, rmon = l.monomer, r.monomer
            lhs = l.site_conditions
            rhs = r.site_conditions
            if lmon.name == rmon.name:
                diffs = [(k, v) for k, v in rhs.items() if lhs[k] != v]
            else:
                diffs = rhs.items()
            key.append((lmon.name, tuple(sorted(lhs.keys())),
                        rmon.name, tuple(sorted(diffs))))
        for r in rh[n:]:
            rhs = r.site_conditions
            key.append((None, tuple(),
                        r.monomer.name, tuple(sorted(rhs.items()))))
        for l in lh[n:]:
            lhs = l.site_conditions
            key.append((l.monomer.name, tuple(sorted(lhs.keys())),
                        None, tuple()))
        grouped[tuple(key)].append(rule)
    return {k: g for k, g in grouped.items() if len(g) > 1}


def _get_lh_rh(rule):
    rexp = rule.rule_expression
    lh = [r.monomer_patterns[0]
          for r in rexp.reactant_pattern.complex_patterns]
    rh = [r.monomer_patterns[0]
          for r in rexp.product_pattern.complex_patterns]
    return lh, rh


def context_elimination(model, sort_key=None):
    '''
    Yields models that are in relation context elimination
    with the model passed in.

    Parameters
    ----------
    model: pysb.Model
    sort_key: function
        Function to be used for sorting the rules.
    '''
    for i, rule in enumerate(model.rules):
        lh, rh = _get_lh_rh(rule)
        n = min(len(lh), len(rh))
        for i in range(n):
            l, r = lh[i], rh[i]
            for site in l.site_conditions.keys():
                if l.monomer == r.monomer and \
                        l.site_conditions[site] == r.site_conditions[site]:
                    mon = l.monomer
                    # site does not change --> is context, try to remove it
                    l_ = mon(**{k: v for k, v in l.site_conditions.items()
                                if k != site})
                    lh_ = reduce(add, lh[:i] + [l_] + lh[i + 1:])
                    r_ = mon(**{k: v for k, v in r.site_conditions.items()
                                if k != site})
                    rh_ = reduce(add, rh[:i] + [r_] + rh[i + 1:])
                    rule_new = inherit_rule(rule, lh_, rh_)
                    m = copy_no_rules(model)
                    # add the rules in preserved order:
                    rules = [rul for rul in model.rules
                             if rul.name != rule_new.name]
                    rules.append(rule_new)
                    rules.sort(key=sort_key)
                    for rul in rules:
                        m.add_component(rul)
                    yield 'ce: %s: %s(%s~%s)' % (rule.name, mon.name, site,
                                                 l.site_conditions[site]), m


def rules_removal(model, sort_key=None):
    '''
    Yields models that are in relation rule elimination
    with the model passed in.

    Parameters
    ----------
    model: pysb.Model
    sort_key: function
        Function to be used for sorting the rules.
    '''
    for rule in model.rules:
        m = copy_no_rules(model)
        for r in model.rules:
            if r != rule:
                m.add_component(r)
        yield 'rem: %s' % rule.name, m


def inherit_rule(old_rule, new_lh, new_rh, name=None):
    '''
    Inherits the rates from the old rule but uses the new
    left- and right-hand sides.

    Parameters
    ----------
    old_rule : pysb.Rule
        The rule to inherit from
    new_lhs : pysb.Expression
        The new left-hand side.
    new_rhs : pysb.Expression
        The new right-hand side.
    name : string or None
        Name for the new rule. If None then the name of the old rule is used.
    '''
    rule = old_rule
    rexp = rule.rule_expression
    rexp_ = new_lh <> new_rh if rexp.is_reversible else new_lh >> new_rh
    return Rule(name or rule.name, rexp_, rule.rate_forward,
                rule.rate_reverse, _export=False)


def get_rule_sort_key(model):
    '''
    Returns a function that can be passed to a sort as a key parameter.
    This sorter tries to keep the rules in the same order as they were in the
    original model.

    Parameters
    ----------
    model: pysb.Model
    '''
    rule_d = {}
    for i, rule in enumerate(model.rules):
        rule_d[rule.name] = i
    return lambda r: rule_d.get(r.name.split('__')[0], len(rule_d))


def parse_reaction_network(rn):
    '''
    Parses a BNGL-generated reaction network into a multiset of edges.
    Useful because the ordering of species and reactions might be different
    althoug the models are equivalent.

    Parameters
    ----------
    rn : string
        A string representing BNG reaction network.
    '''
    # ignore parameters for now
    lines = rn.split('\n')
    species = lines[lines.index('begin species') + 1:
                    lines.index('end species')]
    species_map = {}
    for s in species:
        sid, sp, _ = s.split()
        species_map[sid] = sp
    edges = Counter()
    reactions = lines[lines.index('begin reactions') + 1:
                      lines.index('end reactions')]
    for r in reactions:
        splits = r.split()
        left = ','.join(sorted([species_map[spid]
                                for spid in splits[1].split(',')]))
        right = ','.join(sorted([species_map[spid]
                                 for spid in splits[2].split(',')]))
        rate = splits[3]
        edges.update([(left, right, rate)])
    return edges


def reaction_network(m):
    '''
    Call BNG to construct the reaction network of the model.
    The resulting network is then parsed.

    Parameters
    ----------
    m : pysb.Model
        Model from which the reaction net is to be constructed.
    '''
    return parse_reaction_network(bng.generate_network(m))


def to_bngl(model):
    '''
    Converts the pysb model to BNGL.

    Parameters
    ----------
    model: pysb.Model
    '''
    return bng.BngGenerator(model).get_content()


def rules(model):
    '''
    Return rules of the specified model sorted by their names.
    Useful for canonic model representation.

    Parameters
    ----------
    model: pysb.Model
    '''
    return sorted([r for r in model.rules], key=lambda r: r.name)


def dfs(model, up_to_cee=False):
    '''
    Performs a depth-first search in a space of equivalent models given by
    syntactic operations. Two models are considered equivalent if they
    generate the same reaction network.
    Returns a model that cannot be further simplified by syntactic operations.

    Parameters
    ----------
    model : pysb.Model
        Initial model.
    up_to_cee : bool
        If True, then it performes only bidirectional merging and context
        enumeration elimination.
    '''
    # fast, finds one fix point
    sort_key = get_rule_sort_key(model)
    rn1 = parse_reaction_network(bng.generate_network(model))
    node = model
    while node:
        m1 = node
        node = None
        for edge, m2 in search_model(m1, sort_key=sort_key,
                                     up_to_cee=up_to_cee):
            rn2 = parse_reaction_network(bng.generate_network(m2))
            if rn1 == rn2:
                print edge
                node = m2
                break
    return m1


def copy_no_rules(model):
    '''
    Copies a model without rules.

    Parameters
    ----------
    model : pysb.Model
        Model to copy.
    '''
    m = Model(_export=False)
    for comp in model.all_components():
        if comp.__class__ is not Rule:
            m.add_component(comp)
    for ini in model.initial_conditions:
        m.initial(*ini)
    return m


from pysb.bng import BngGenerator, _get_bng_path, _parse_bng_outfile
from pysb.bng import GenerateNetworkError
import os
import random
import subprocess


def bng_simulate(model, times, method='ode', output_dir='/tmp', cleanup=True):
    """
    Simulate a model with BNG's simulator and return the trajectories.
    Adapted from pysb.bng.run_ssa.

    Parameters
    ----------
    model : pysb.Model or string
        Model to simulate. Can be either a pysb.Model or a string representing
        BNGL model.
    times: list of floats
        Sample times.
    method: string
        'ode' or 'ssa'
    output_dir : string, optional
        Location for temporary files generated by BNG. Defaults to '/tmp'.
    cleanup : bool, optional
        If True (default), delete the temporary files after the simulation is
        finished. If False, leave them in place (in `output_dir`). Useful for
        debugging.

    """

    times = list(times)
    run_ssa_code = """
    begin actions
    generate_network({overwrite=>1});
    simulate_%s({sample_times=>%s});\n
    end actions
    """ % (method, times)

    if not isinstance(model, str):
        model = BngGenerator(model).get_content()
    bng_filename = '%d_%d_temp.bngl' % (
        os.getpid(), random.randint(0, 10000))
    gdat_filename = bng_filename.replace('.bngl', '.gdat')
    cdat_filename = bng_filename.replace('.bngl', '.cdat')
    net_filename = bng_filename.replace('.bngl', '.net')

    try:
        working_dir = os.getcwd()
        os.chdir(output_dir)
        bng_file = open(bng_filename, 'w')
        bng_file.write(model)
        bng_file.write(run_ssa_code)
        bng_file.close()
        p = subprocess.Popen(['perl', _get_bng_path(), bng_filename],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (p_out, p_err) = p.communicate()
        if p.returncode:
            raise GenerateNetworkError(p_out.rstrip("at line") + "\n" +
                                       p_err.rstrip())

        output_arr = _parse_bng_outfile(gdat_filename)
    finally:
        if cleanup:
            for filename in [bng_filename, gdat_filename,
                             cdat_filename, net_filename]:
                if os.access(filename, os.F_OK):
                    os.unlink(filename)
        os.chdir(working_dir)
    return output_arr
