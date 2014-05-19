from pysb import Model, Rule, bng
from collections import defaultdict, Counter


def search_model(model, sort_key=None):
    '''
    Yields all descsendants of a model. They are created by applying syntactic
    operations.
    '''
    for m in merge_bidirectional(model, sort_key):
        yield m
    for m in context_enumeration_elimination(model, sort_key):
        yield m
    for m in context_elimination(model, sort_key):
        yield m
    for m in unreachable_rules_removal(model):
        yield m


def merge_bidirectional(model, sort_key=None):
    lh_map = defaultdict(list)
    rh_map = defaultdict(list)
    for rule in model.rules:
        if not rule.rule_expression.is_reversible:
            lh, rh = get_lh_rh(rule)
            lh_map[(repr(lh), repr(rh))].append(rule)
            rh_map[(repr(rh), repr(lh))].append(rule)
    for rep in set(lh_map.keys()) & set(rh_map.keys()):
        for fwd in lh_map[rep]:
            for bwd in rh_map[rep]:
                m = copy_no_rules(model)
                rules = [r for r in model.rules if r not in [fwd, bwd]]
                lh, rh = get_lh_rh(fwd)
                new_rule = Rule('__'.join([fwd.name, bwd.name]), lh <> rh,
                                fwd.rate_forward, bwd.rate_forward,
                                _export=False)
                rules.append(new_rule)
                rules.sort(key=sort_key)
                for r in rules:
                    m.add_component(r)
                yield 'bi: %s' % new_rule.name, m


def context_enumeration_elimination(model, sort_key=None):
    grouped = _group_rules_by_context_and_changes(model.rules)
    for (_, _, diffs), group in grouped.items():
        mon = group[0].reactant_pattern.complex_patterns[0].\
            monomer_patterns[0].monomer
        diffs = set(zip(*diffs)[0])
        sites = group[0].reactant_pattern.complex_patterns[0].\
            monomer_patterns[0].site_conditions.keys()
        sites = [s for s in sites if s not in diffs]
        unchanged = []
        for rule in group:
            lh, _ = get_lh_rh(rule)
            lhs = lh.site_conditions
            unchanged.append(({k: v for k, v in lhs.items() if k not in diffs},
                              rule))
        for site in sites:
            # local map states of the `site` into sets of occuring contexts
            local = defaultdict(set)
            # conmap maps contexts to lists of rules that contain that context
            conmap = defaultdict(list)
            excludes = set(list(diffs) + [site])
            for rule in group:
                lhs = get_lh_rh(rule)[0].site_conditions
                context = tuple(sorted([(k, v) for k, v in lhs.items()
                                        if k not in excludes]))
                local[lhs[site]].add(context)
                conmap[context].append(rule)
            ref = local[local.keys()[0]]
            possible_states = mon.site_states[site]
            if all(local[st] == ref for st in possible_states):
                #  it is enumerating, we can reduce :)
                m = copy_no_rules(model)
                rules = [r for r in model.rules if r not in group]
                new_rules = []
                for context in ref:
                    name = '__'.join([r.name for r in conmap[context]])
                    rule = conmap[context][0]
                    rexp = rule.rule_expression
                    # now set the site state to None
                    lh = rexp.reactant_pattern.complex_patterns[0].copy()
                    rh = rexp.product_pattern.complex_patterns[0].copy()
                    lh.monomer_patterns[0].site_conditions[site] = None
                    rh.monomer_patterns[0].site_conditions[site] = None
                    rule_new = inherit_rule(rule, lh, rh, name=name)
                    new_rules.append(rule_new)
                rules.extend(new_rules)
                rules.sort(key=sort_key)
                for r in rules:
                    m.add_component(r)
                yield 'CEE: %s %s; %s' % (
                    mon.name, site,
                    ' '.join(map(lambda r: r.name, new_rules))), m


def _group_rules_by_context_and_changes(rules):
    grouped = defaultdict(list)
    for rule in rules:
        lh, rh = get_lh_rh(rule)
        mon = lh.monomer
        lhs = lh.site_conditions
        rhs = rh.site_conditions
        diffs = [(k, v) for k, v in rhs.items() if lhs[k] != v]
        key = (mon.name, tuple(sorted(lhs.keys())), tuple(sorted(diffs)))
        grouped[key].append(rule)
    return {k: g for k, g in grouped.items() if len(g) > 1}


def get_lh_rh(rule):
    rexp = rule.rule_expression
    react = rexp.reactant_pattern.complex_patterns
    if len(react) != 1:
        raise ValueError("Can't handle this yet - too many reactants")
    prod = rexp.product_pattern.complex_patterns
    if len(prod) != 1:
        raise ValueError("Can't handle this yet - too many products")
    lh, rh = react[0].monomer_patterns[0], prod[0].monomer_patterns[0]
    if lh.monomer != rh.monomer:
        raise ValueError("LHS and RHS do not match")
    if set(lh.site_conditions.keys()) != set(rh.site_conditions.keys()):
        raise ValueError("LHS and RHS sites do not match")
    return lh, rh


def context_elimination(model, sort_key=None):
    for i, rule in enumerate(model.rules):
        lh, rh = get_lh_rh(rule)
        mon = lh.monomer
        for site in lh.site_conditions.keys():
            if lh.site_conditions[site] == rh.site_conditions[site]:
                # site does not change --> is context, try to remove it
                lh_ = mon(**{k: v for k, v in lh.site_conditions.items()
                             if k != site})
                rh_ = mon(**{k: v for k, v in rh.site_conditions.items()
                             if k != site})
                rule_new = inherit_rule(rule, lh_, rh_)
                m = copy_no_rules(model)
                # add the rules in preserved order:
                rules = [r for r in model.rules if r.name != rule_new.name]
                rules.append(rule_new)
                rules.sort(key=sort_key)
                for r in rules:
                    m.add_component(r)
                yield 'ce: %s: %s(%s~%s)' % (rule.name, mon.name, site,
                                             lh.site_conditions[site]), m


def unreachable_rules_removal(model, sort_key=None):
    for rule in model.rules:
        m = copy_no_rules(model)
        for r in model.rules:
            if r != rule:
                m.add_component(r)
        yield 'rem: %s' % rule.name, m


def inherit_rule(old_rule, new_lh, new_rh, name=None):
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
        left = species_map[splits[1]]
        right = species_map[splits[2]]
        rate = splits[3]
        edges.update([(left, right, rate)])
    return edges


def reaction_network(m):
    return parse_reaction_network(bng.generate_network(m))


def to_bngl(model):
    return bng.BngGenerator(model).get_content()


def rules(model):
    '''
    Return rules of the specified model sorted by their names.
    Useful for canonic model representation.
    '''
    return sorted([r for r in model.rules], key=lambda r: r.name)


def dfs(model):
    '''
    Performs a depth-first search in a space of equivalent models given by
    syntactic operations. Two models are considered equivalent if they
    generate the same reaction network.
    Returns a model that cannot be further simplified by syntactic operations.
    '''
    # fast, finds one fix point
    sort_key = get_rule_sort_key(model)
    rn1 = parse_reaction_network(bng.generate_network(model))
    node = model
    while node:
        m1 = node
        node = None
        for edge, m2 in search_model(m1, sort_key=sort_key):
            rn2 = parse_reaction_network(bng.generate_network(m2))
            if rn1 == rn2:
                print edge
                node = m2
                break
    return m1


def copy_no_rules(model):
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
    model : Model
        Model to simulate.
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

    gen = BngGenerator(model)
    bng_filename = '%s_%d_%d_temp.bngl' % (
        model.name, os.getpid(), random.randint(0, 10000))
    gdat_filename = bng_filename.replace('.bngl', '.gdat')
    cdat_filename = bng_filename.replace('.bngl', '.cdat')
    net_filename = bng_filename.replace('.bngl', '.net')

    try:
        working_dir = os.getcwd()
        os.chdir(output_dir)
        bng_file = open(bng_filename, 'w')
        bng_file.write(gen.get_content())
        bng_file.write(run_ssa_code)
        bng_file.close()
        p = subprocess.Popen(['perl', _get_bng_path(), bng_filename],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (p_out, p_err) = p.communicate()
        if p.returncode:
            raise GenerateNetworkError(p_out.rstrip("at line") + "\n" +
                                       p_err.rstrip())

        output_arr = _parse_bng_outfile(gdat_filename)
        #ssa_file = open(ssa_filename, 'r')
        #output.write(ssa_file.read())
        #net_file.close()
        #if append_stdout:
        #    output.write("#\n# BioNetGen execution log follows\n# ==========")
        #    output.write(re.sub(r'(^|\n)', r'\n# ', p_out))
    finally:
        if cleanup:
            for filename in [bng_filename, gdat_filename,
                             cdat_filename, net_filename]:
                if os.access(filename, os.F_OK):
                    os.unlink(filename)
        os.chdir(working_dir)
    return output_arr
