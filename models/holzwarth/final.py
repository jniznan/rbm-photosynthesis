from pysb import *

Model()

param_dict = {
    'deex_flu': 670e5,  # de-excitation (fluorescence)
    'deex_heat': 2e5,  # de-excitation (heat)
    'ch_re': 50e10,  # Charge recombination
    'ch_se': 50e10,  # Primary charge separation
    'stab_gen': 25e10,  # ChlD1+ reduction by P680
    'stab_deg': 6.7e10,  # Reverse electron flow from ChlD1 to P680+
    'qa_red': 0.48e10,  # Qa reduction by PhD1-
    'qa_oxi': 0.24e10,  # Reverse electron flow from Qa- to PhD1
    'extrans': 1.92e10,  # Energy transfer between A and ChlD1  k1
    'extrans_re': 2.5e10,  # Energy transfer between A and ChlD1  k2
}

for p, v in param_dict.items():
    Parameter(p, v)

Monomer('PSII', ['P680', 'Qa', 'PhD1', 'ac', 'ChlD1'],
        {'P680': ['n', 'p'], 'Qa': ['n', 'm'], 'PhD1': ['n', 'm'],
         'ac': ['n', 'exc'], 'ChlD1': ['n', 'p', 'exc']})

# rules
Rule('deexcitation_antenna_fluorescence',
     PSII(ac='exc') >> PSII(ac='n'),
     deex_flu)
Rule('deexcitation_antenna_heat',
     PSII(ac='exc') >> PSII(ac='n'),
     deex_heat)
Rule('deexcitation_ChlD1_fluorescence',
     PSII(ChlD1='exc') >> PSII(ChlD1='n'),
     deex_flu)
Rule('deexcitation_ChlD1_heat',
     PSII(ChlD1='exc') >> PSII(ChlD1='n'),
     deex_heat)
Rule('excitation_transfer',
     PSII(ac='exc', ChlD1='n') <> PSII(ac='n', ChlD1='exc'),
     extrans, extrans_re)
Rule('primary_charge_separation_recombination',
     PSII(PhD1='n', ChlD1='exc') <> PSII(PhD1='m', ChlD1='p'),
     ch_se, ch_re)
Rule('stable_pair_generation_degeneration',
     PSII(P680='n', PhD1='m', ChlD1='p') <>
     PSII(P680='p', PhD1='m', ChlD1='n'),
     stab_gen, stab_deg)
Rule('quinone_qa_reduction_oxidation',
     PSII(P680='p', Qa='n', PhD1='m') <>
     PSII(P680='p', Qa='m', PhD1='n'),
     qa_red, qa_oxi)

# initial conditions
Parameter('init', 1)
Initial(PSII(P680='n', Qa='n', PhD1='n', ac='exc', ChlD1='n'), init)

# observables
Observable('Ant_exc', PSII(P680='n', Qa='n', PhD1='n', ac='exc', ChlD1='n'))
Observable('ChlD1_exc', PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='exc'))
