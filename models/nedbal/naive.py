from pysb import *

Model()

# Since this model does nt have the rates available, we create here
# a dummy rate parameter in order to be able to express the model in
# rule-based form.
Parameter('k', 1.)

Monomer('PSII', ['Ant', 'Yz', 'P680', 'ChlD', 'Pheo', 'Qa', 'S'],
        {'Ant': ['n', 'exc'],
         'Yz': ['n', 'p'],
         'P680': ['n', 'p'],
         'ChlD': ['n', 'exc', 'p'],
         'Pheo': ['n', 'm'],
         'Qa': ['n', 'm'],
         'S': ['0', '1', '2', '3', '4']})

# rules
Rule('excitation_antenna',
     PSII(Ant='n') >> PSII(Ant='exc'), k)
Rule('excitation_primary_chlorophyll_donor',
     PSII(ChlD='n') >> PSII(ChlD='exc'), k)
Rule('energy_transfer',
     PSII(Ant='exc', ChlD='n') <> PSII(Ant='n', ChlD='exc'), k, k)
Rule('deexcitation_antenna_fluorescence',
     PSII(Ant='exc') >> PSII(Ant='n'), k)
Rule('deexcitation_chlor_donor_fluorescence',
     PSII(ChlD='exc') >> PSII(ChlD='n'), k)
Rule('primary_charge_separation',
     PSII(Yz='n', P680='n', ChlD='exc', Pheo='n') >>
     PSII(Yz='n', P680='n', ChlD='p', Pheo='m'), k)
Rule('stable_radical_pair_generation',
     PSII(Yz='n', P680='n', ChlD='p', Pheo='m') >>
     PSII(Yz='n', P680='p', ChlD='n', Pheo='m'), k)
Rule('quinone_qa_reduction_preserving_P680_oxidized',
     PSII(Yz='n', P680='p', ChlD='n', Pheo='m', Qa='n') >>
     PSII(Yz='n', P680='p', ChlD='n', Pheo='n', Qa='m'), k)
Rule('tyrosine_oxidation',
     PSII(Yz='n', P680='p', ChlD='n', Pheo='n', Qa='m') >>
     PSII(Yz='p', P680='n', ChlD='n', Pheo='n', Qa='m'), k)
Rule('tyrosine_oxidation_preserving_Pheo_reduced',
     PSII(Yz='n', P680='p', ChlD='n', Pheo='m') >>
     PSII(Yz='p', P680='n', ChlD='n', Pheo='m'), k)
Rule('quinone_qa_reduction_preserving_Yz_oxidized',
     PSII(Yz='p', P680='n', ChlD='n', Pheo='m', Qa='n') >>
     PSII(Yz='p', P680='n', ChlD='n', Pheo='n', Qa='m'), k)
Rule('quinone_qa_reduction',
     PSII(Yz='n', P680='n', ChlD='n', Pheo='m', Qa='n') >>
     PSII(Yz='n', P680='n', ChlD='n', Pheo='n', Qa='m'), k)
for s in range(4):
    Rule('s_state_transition_%s' % s,
         PSII(Yz='p', P680='n', ChlD='n', Pheo='n', S=str(s)) >>
         PSII(Yz='n', P680='n', ChlD='n', Pheo='n', S=str(s + 1)), k)
for s in range(4):
    Rule('s_state_transition_%s_with_Pheo_reduced' % s,
         PSII(Yz='p', P680='n', ChlD='n', Pheo='m', S=str(s)) >>
         PSII(Yz='n', P680='n', ChlD='n', Pheo='m', S=str(s + 1)), k)
Rule('S4_to_S0_transition',
     PSII(S='4') >> PSII(S='0'), k)  # TODO include oxygen release ??


# initial conditions
Parameter('init', 1)
Initial(PSII(Ant='n', Yz='n', P680='n', ChlD='n', Pheo='n', Qa='n', S='0'),
        init)
