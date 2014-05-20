
from pysb import *

Model()

param_dict = {
    'kHCL': 7,
    'k01': 20000,
    'k12': 10000,
    'k23': 3330,
    'k30': 1000,
    'ka': 0.01,
    'kAB1': 3500,
    'kAB2': 1750,
    'kbHC': 10,
    'kBA1': 175,
    'kBA2': 35,
    'kbB': 250,
    'kbct': 100,
    'kbF': 10,
    'kbFd': 5,
    'kbLF': 10,
    'kbR': 10,
    'kbX': 10,
    'kfHC': 100,
    'kfB': 250,
    'kfct': 100,
    'kfF': 100,
    'kfFd': 5,
    'kfLF': 100,
    'kFNR': 220,
    'kfR': 100,
    'kfX': 100,
    'kLHC': 2300,
    'kL1': 200,
    'kL2': 4000,
}


for p, v in param_dict.items():
    Parameter(p, v)


Monomer('PSII', ['P680', 'Qa', 'Qb'],
        {'P680': ['n', 'p'], 'Qa': ['n', 'm'], 'Qb': ['n', 'm', '2m']})
Monomer('PQ', [], {})
Monomer('PQH', [], {})
Monomer('S', ['x'], {'x': ['0', '1', '2', '3']})
Monomer('CytB6F', ['bL', 'bHc', 'f'],
        {'bL': ['n', 'm'], 'bHc': ['n', 'm', '2m'], 'f': ['n', 'm']})
Monomer('Fd', ['x'], {'x': ['n', 'm']})
Monomer('Pc', ['x'], {'x': ['n', 'p']})
Monomer('PSI', ['P700', 'Fb'], {'P700': ['n', 'p'], 'Fb': ['n', 'm']})
Monomer('FNR', ['x'], {'x': ['i', 'a', 'am', 'a2m']})


Rule('electron_transport_from_Fd_to_FNR_1', FNR(x='a') + Fd(x='m') <> FNR(x='am') + Fd(x='n'), kfFd, kbFd),
Rule('electron_transport_from_Fd_to_FNR_2', FNR(x='am') + Fd(x='m') <> FNR(x='a2m') + Fd(x='n'), kfFd, kbFd),
Rule('activation_of_FNR', FNR(x='i') >> FNR(x='a'), ka),
Rule('turnover_of_FNR', FNR(x='a2m') >> FNR(x='a'), kFNR),
Rule('charge_separation_in_PSI', PSI(P700='n', Fb='n') >> PSI(P700='p', Fb='m'), kL1),
Rule('electron_transport_from_photosystemI_to_Fd___1__2', Fd(x='n') + PSI(Fb='m') <> Fd(x='m') + PSI(Fb='n'), kfX, kbX),
Rule('electron_transports_from_qa_to_qb___1__2', PSII(Qa='m', Qb='n') <> PSII(Qa='n', Qb='m'), kAB1, kBA1),
Rule('electron_transports_from_qa_to_reduced_qb___1__2', PSII(Qa='m', Qb='m') <> PSII(Qa='n', Qb='2m'), kAB2, kBA2),
Rule('electron_transport_inside_cytochrome___3__4', CytB6F(bL='m', bHc='m') <> CytB6F(bL='n', bHc='2m'), kLHC, kHCL),
Rule('electron_transport_inside_cytochrome___1__2', CytB6F(bL='m', bHc='n') <> CytB6F(bL='n', bHc='m'), kLHC, kHCL),
Rule('electron_transport_from_reduced_PQ_to_cytochrome___1__2__3', PQH() + CytB6F(bL='n', f='n') <> PQ() + CytB6F(bL='m', f='m'), kfLF, kbLF),
Rule('charge_separation_in_PSII___1__2__3', PSII(P680='n', Qa='n') >> PSII(P680='p', Qa='m'), kL2),
Rule('electron_donation_to_P700_by_Pc___1__2', Pc(x='n') + PSI(P700='p') <> Pc(x='p') + PSI(P700='n'), kfR, kbR),
Rule('electron_donation___1_2__2_2__3_2__4_2__5_2__6_2', S(x='1') + PSII(P680='p') >> S(x='2') + PSII(P680='n'), k12),
Rule('electron_transport_from_cytochrome_to_oxidized_PQ___1__2__3__4', PQ() + CytB6F(bHc='2m') <> PQH() + CytB6F(bHc='n'), kfHC, kbHC),
Rule('electron_transport_from_reduced_Fd_to_cytochrome___1__2__3__4', Fd(x='m') + CytB6F(bHc='n') <> Fd(x='n') + CytB6F(bHc='m'), kfct, kbct),
Rule('electron_donation___1_1__2_1__3_1__4_1__5_1__6_1', S(x='0') + PSII(P680='p') >> S(x='1') + PSII(P680='n'), k01),
Rule('electron_donation___1_3__2_3__3_3__4_3__5_3__6_3', S(x='2') + PSII(P680='p') >> S(x='3') + PSII(P680='n'), k23),
Rule('electron_transport_from_reduced_Fd_to_cytochrome___5__6__7__8', Fd(x='m') + CytB6F(bHc='m') <> Fd(x='n') + CytB6F(bHc='2m'), kfct, kbct),
Rule('Qb_PQ_exchange___1__2__3__4', PQ() + PSII(Qb='2m') <> PQH() + PSII(Qb='n'), kfB, kbB),
Rule('electron_donation___1_0__2_0__3_0__4_0__5_0__6_0', S(x='3') + PSII(P680='p') >> S(x='0') + PSII(P680='n'), k30),
Rule('electron_transport_from_cytochrome_to_oxidized_Pc___1__3__5__2__4__6', CytB6F(f='m') + Pc(x='p') <> CytB6F(f='n') + Pc(x='n'), kfF, kbF),


# initial conditions
#bL/bH.c/f   1
Parameter('init_cyt', 1)
Initial(CytB6F(bL='n', bHc='n', f='n'), init_cyt)
#Fd  3
Parameter('init_Fd', 3)
Initial(Fd(x='n'), init_Fd)
#FNRi    3
Parameter('init_FNR', 3)
Initial(FNR(x='i'), init_FNR)
#P680/Qa/Qb  1
Parameter('init_ps2', 1)
Initial(PSII(P680='n', Qa='n', Qb='n'), init_ps2)
#P700/Fb 1
Parameter('init_ps1', 1)
Initial(PSI(P700='n', Fb='n'), init_ps1)
#Pc  3
Parameter('init_Pc', 3)
Initial(Pc(x='n'), init_Pc)
#PQ  2.5
Parameter('init_PQ', 2.5)
Initial(PQ(), init_PQ)
#PQH 2.5
Parameter('init_PQH', 2.5)
Initial(PQH(), init_PQH)
#S0  0.25
Parameter('init_s0', 0.25)
Initial(S(x='0'), init_s0)
#S1  0.75
Parameter('init_s1', 0.75)
Initial(S(x='1'), init_s1)


# observables
Observable('Qa_m', PSII(Qa='m'))
Observable('PQ_n', PQ())
Observable('Pc_p', Pc(x='p'))
Observable('P700_p', PSI(P700='p'))
