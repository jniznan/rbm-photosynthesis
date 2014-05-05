from pysb import Model, Rule, bng
from collections import defaultdict, Counter


def search_model(model):
    '''
    Yields all descsendants of a model. They are created by applying syntactic
    operations.
    '''
    for m in context_enumeration_elimination(model):
        yield m
    for m in context_elimination(model):
        yield m


def merge_bidirectional(model):
    pass


def context_enumeration_elimination(model):
    sort_key = get_rule_sort_key(model)
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


def context_elimination(model):
    for i, rule in enumerate(model.rules):
        lh, rh = get_lh_rh(rule)
        mon = lh.monomer
        sort_key = get_rule_sort_key(model)
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


def rules(model):
    '''
    Return rules of the specified model sorted by their names.
    Useful for canonic model representation.
    '''
    return sorted([r for r in model.rules], key=lambda r: r.name)


def bfs(model):
    '''
    Performs a breadth-first search in a space of equivalent models given by
    syntactic operations. Two models are considered equivalent if they
    generate the same reaction network.
    '''
    # WARNING: takes looong time, finds all fix points
    rn1 = parse_reaction_network(bng.generate_network(model))

    searched = set()
    fringe = [model]
    edges = []

    while fringe:
        m1 = fringe.pop(0)
        mr1 = repr(rules(m1))
        searched.add(mr1)
        for edge, m2 in search_model(m1):
            rn2 = parse_reaction_network(bng.generate_network(m2))
            if rn1 == rn2:
                mr2 = repr(rules(m2))
                edges.append((mr1, edge, mr2))
                if mr2 not in searched:
                    fringe.append(m2)


def dfs(model):
    '''
    Performs a depth-first search in a space of equivalent models given by
    syntactic operations. Two models are considered equivalent if they
    generate the same reaction network.
    Returns a model that cannot be further simplified by syntactic operations.
    '''
    # fast, finds one fix point
    rn1 = parse_reaction_network(bng.generate_network(model))
    node = model
    while node:
        m1 = node
        node = None
        for edge, m2 in search_model(m1):
            rn2 = parse_reaction_network(bng.generate_network(m2))
            if rn1 == rn2:
                node = m2
                break
    return m1
    #return remove_duplicate_rules(m1)


def remove_duplicate_rules(model):
    m = copy_no_rules(model)
    rule_map = defaultdict(list)
    for rule in model.rules:
        rid = repr(rule.rule_expression) + repr(rule.rate_forward) +\
            repr(rule.rate_reverse)
        rule_map[rid].append(rule)
    # leave only one rule in each group:
    rule_map = {k: min(rules, key=lambda r: (len(r.name), r.name))
                for k, rules in rule_map.items()}
    for r in sorted(rule_map.values(), key=get_rule_sort_key(model)):
        m.add_component(r)
    return m


def copy_no_rules(model):
    m = Model(_export=False)
    for comp in model.all_components():
        if comp.__class__ is not Rule:
            m.add_component(comp)
    for ini in model.initial_conditions:
        m.initial(*ini)
    return m
