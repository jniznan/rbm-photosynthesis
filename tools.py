from pysb import Model, Rule, bng


def search_model(model):
    for i, rule in enumerate(model.rules):
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
        mon = lh.monomer
        if set(lh.site_conditions.keys()) != set(rh.site_conditions.keys()):
            raise ValueError("LHS and RHS site conditions do not match")
        for site in lh.site_conditions.keys():
            if lh.site_conditions[site] == rh.site_conditions[site]:
                # site does not change --> is context, try to remove it
                lh_ = mon(**{k: v for k, v in lh.site_conditions.items()
                             if k != site})
                rh_ = mon(**{k: v for k, v in rh.site_conditions.items()
                             if k != site})
                rexp_ = lh_ <> rh_ if rexp.is_reversible else lh_ >> rh_
                m = Model(_export=False)
                r = Rule(rule.name, rexp_,
                         rule.rate_forward, rule.rate_reverse, _export=False)
                for comp in model.all_components():
                    if comp != rule:
                        m.add_component(comp)
                for ini in model.initial_conditions:
                    m.initial(*ini)
                m.add_component(r)
                yield '%s: %s(%s~%s)' % (rule.name, mon.name, site,
                                         lh.site_conditions[site]), m


def parse_reaction_network(rn):
    # ignore parameters for now
    lines = rn.split('\n')
    species = lines[lines.index('begin species') + 1:
                    lines.index('end species')]
    species_map = {}
    for s in species:
        sid, sp, _ = s.split()
        species_map[sid] = sp
    edges = set()
    reactions = lines[lines.index('begin reactions') + 1:
                      lines.index('end reactions')]
    for r in reactions:
        splits = r.split()
        left = species_map[splits[1]]
        right = species_map[splits[2]]
        edges.add((left, right))
    return edges


def rules(model):
    return sorted([r for r in model.rules], key=lambda r: r.name)


def bfs(model):
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
