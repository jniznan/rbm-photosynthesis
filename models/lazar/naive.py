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
    'kf(HC)': 100,
    'kfB': 250,
    'kfct': 100,
    'kfF': 100,
    'kfFd': 5,
    'kfLF': 100,
    'kFNR': 220,
    'kfR': 100,
    'kfX': 100,
    'kL(HC)': 2300,
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


i = 0
for a in ['n', 'm', '2m']:
    for b in ['n', 'm']:
        i += 1
        Rule('electron_transport_from_cytochrome_to_oxidized_Pc_%s' % i,
             CytB6F(bL=b, bHc=a, f='m') + Pc(x='p') <>
             CytB6F(bL=b, bHc=a, f='n') + Pc(x='n'),
             kfF, kbF)

Rule('electron_transport_from_reduced_Fd_to_cytochrome_1',
     Fd(x='m') + CytB6F(bL='n', bHc='n', f='n') <>
     Fd(x='n') + CytB6F(bL='n', bHc='m', f='n'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_2',
     Fd(x='m') + CytB6F(bL='m', bHc='n', f='n') <>
     Fd(x='n') + CytB6F(bL='m', bHc='m', f='n'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_3',
     Fd(x='m') + CytB6F(bL='n', bHc='n', f='m') <>
     Fd(x='n') + CytB6F(bL='n', bHc='m', f='m'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_4',
     Fd(x='m') + CytB6F(bL='m', bHc='n', f='m') <>
     Fd(x='n') + CytB6F(bL='m', bHc='m', f='m'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_5',
     Fd(x='m') + CytB6F(bL='n', bHc='m', f='n') <>
     Fd(x='n') + CytB6F(bL='n', bHc='2m', f='n'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_6',
     Fd(x='m') + CytB6F(bL='m', bHc='m', f='n') <>
     Fd(x='n') + CytB6F(bL='m', bHc='2m', f='n'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_7',
     Fd(x='m') + CytB6F(bL='n', bHc='m', f='m') <>
     Fd(x='n') + CytB6F(bL='n', bHc='2m', f='m'),
     kfct, kbct)
Rule('electron_transport_from_reduced_Fd_to_cytochrome_8',
     Fd(x='m') + CytB6F(bL='m', bHc='m', f='m') <>
     Fd(x='n') + CytB6F(bL='m', bHc='2m', f='m'),
     kfct, kbct)

Rule('electron_transport_inside_cytochrome_1',
     CytB6F(bL='m', bHc='n', f='n') <> CytB6F(bL='n', bHc='m', f='n'),
     kLHC, kHCL)
Rule('electron_transport_inside_cytochrome_2',
     CytB6F(bL='m', bHc='n', f='m') <> CytB6F(bL='n', bHc='m', f='m'),
     kLHC, kHCL)
Rule('electron_transport_inside_cytochrome_3',
     CytB6F(bL='m', bHc='m', f='n') <> CytB6F(bL='n', bHc='2m', f='n'),
     kLHC, kHCL)
Rule('electron_transport_inside_cytochrome_4',
     CytB6F(bL='m', bHc='m', f='m') <> CytB6F(bL='n', bHc='2m', f='m'),
     kLHC, kHCL)

Rule('electron_transport_from_reduced_PQ_to_cytochrome_1',
     PQH() + CytB6F(bL='n', bHc='n', f='n') <>
     PQ() + CytB6F(bL='m', bHc='n', f='m'),
     kfLF, kbLF)
Rule('electron_transport_from_reduced_PQ_to_cytochrome_2',
     PQH() + CytB6F(bL='n', bHc='m', f='n') <>
     PQ() + CytB6F(bL='m', bHc='m', f='m'),
     kfLF, kbLF)
Rule('electron_transport_from_reduced_PQ_to_cytochrome_3',
     PQH() + CytB6F(bL='n', bHc='2m', f='n') <>
     PQ() + CytB6F(bL='m', bHc='2m', f='m'),
     kfLF, kbLF)

Rule('electron_transport_from_cytochrome_to_oxidized_PQ_1',
     PQ() + CytB6F(bL='n', bHc='2m', f='n') <>
     PQH() + CytB6F(bL='n', bHc='n', f='n'),
     kfHC, kbHC)
Rule('electron_transport_from_cytochrome_to_oxidized_PQ_2',
     PQ() + CytB6F(bL='m', bHc='2m', f='n') <>
     PQH() + CytB6F(bL='m', bHc='n', f='n'),
     kfHC, kbHC)
Rule('electron_transport_from_cytochrome_to_oxidized_PQ_4',
     PQ() + CytB6F(bL='m', bHc='2m', f='m') <>
     PQH() + CytB6F(bL='m', bHc='n', f='m'),
     kfHC, kbHC)
Rule('electron_transport_from_cytochrome_to_oxidized_PQ_3',
     PQ() + CytB6F(bL='n', bHc='2m', f='m') <>
     PQH() + CytB6F(bL='n', bHc='n', f='m'),
     kfHC, kbHC)

Rule('electron_transport_from_Fd_to_FNR_1',
     FNR(x='a') + Fd(x='m') <> FNR(x='am') + Fd(x='n'),
     kfFd, kbFd)
Rule('electron_transport_from_Fd_to_FNR_2',
     FNR(x='am') + Fd(x='m') <> FNR(x='a2m') + Fd(x='n'),
     kfFd, kbFd)

Rule('electron_transport_from_photosystemI_to_Fd_1',
     Fd(x='n') + PSI(P700='n', Fb='m') <> Fd(x='m') + PSI(P700='n', Fb='n'),
     kfX, kbX)
Rule('electron_transport_from_photosystemI_to_Fd_2',
     Fd(x='n') + PSI(P700='p', Fb='m') <> Fd(x='m') + PSI(P700='p', Fb='n'),
     kfX, kbX)

Rule('activation_of_FNR',
     FNR(x='i') >> FNR(x='a'), ka)

Rule('turnover_of_FNR',
     FNR(x='a2m') >> FNR(x='a'), kFNR)

for a, b, rate in zip(['0', '1', '2', '3'], ['1', '2', '3', '0'],
                      [k01, k12, k23, k30]):
    i = 0
    for qb in ['n', 'm', '2m']:
        for qa in ['n', 'm']:
            i += 1
            Rule('electron_donation_%s_%s' % (i, b),
                 S(x=a) + PSII(P680='p', Qa=qa, Qb=qb) >>
                 S(x=b) + PSII(P680='n', Qa=qa, Qb=qb),
                 rate)


Rule('electron_transports_from_qa_to_qb_1',
     PSII(P680='p', Qa='m', Qb='n') <> PSII(P680='p', Qa='n', Qb='m'),
     kAB1, kBA1)
Rule('electron_transports_from_qa_to_qb_2',
     PSII(P680='n', Qa='m', Qb='n') <> PSII(P680='n', Qa='n', Qb='m'),
     kAB1, kBA1)

i = 0
for p680 in ['p', 'n']:
    i += 1
    Rule('electron_transports_from_qa_to_reduced_qb_%s' % i,
         PSII(P680=p680, Qa='m', Qb='m') <> PSII(P680=p680, Qa='n', Qb='2m'),
         kAB2, kBA2)

i = 0
for qa in ['n', 'm']:
    for p680 in ['p', 'n']:
        i += 1
        Rule('Qb_PQ_exchange_%s' % i,
             PQ() + PSII(P680=p680, Qa=qa, Qb='2m') <>
             PQH() + PSII(P680=p680, Qa=qa, Qb='n'),
             kfB, kbB)

for i, fb in enumerate(['n', 'm']):
    Rule('electron_donation_to_P700_by_Pc_%s' % (i + 1),
         Pc(x='n') + PSI(P700='p', Fb=fb) <>
         Pc(x='p') + PSI(P700='n', Fb=fb),
         kfR, kbR)

Rule('charge_separation_in_PSI',
     PSI(P700='n', Fb='n') > PSI(P700='p', Fb='m'),
     kL1)

for i, qb in enumerate(['n', 'm', '2m']):
    Rule('charge_separation_in_PSII_%s' % (i + 1),
         PSII(P680='n', Qa='n', Qb=qb) <>
         PSII(P680='p', Qa='m', Qb=qb),
         kL2)

#+Electron donation 1_5
#Reaction rate: k01*S0 * P680+/Qa-/N
#Product: S1, P680/Qa-/N
#
#+Electron donation 4_5
#Reaction rate: k1*S2 * P680+/Qa-/N
#Product: S0, P680/Qa-/N
#Kinetic rate constant   Value
#k1  1000


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
Initial(PSI(P700='n', Fb='n'))
#Pc  3
Parameter('init_Pc', 3)
Initial(Pc(x='n'))
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
Initial(S(x='1'))


Assigned quantities
+Fq(t)
Initial expression: ((1-0.55)*("P680+/Qa-/Qb"+"P680+/Qa-/Qb-"+"P680+/Qa-/Qb2-"+"P680/Qa-/Qb2-"+"P680/Qa-/Qb"+"P680/Qa-/Qb-")/(1-0.55*("P680+/Qa-/Qb"+"P680+/Qa-/Qb-"+"P680+/Qa-/Qb2-"+"P680/Qa-/Qb2-"+"P680/Qa-/Qb"+"P680/Qa-/Qb-")))/(1+((1/45+(("P680+/Qa-/Qb"+"P680+/Qa-/Qb-"+"P680+/Qa-/Qb2-"+"P680/Qa-/Qb2-"+"P680/Qa-/Qb"+"P680/Qa-/Qb-")*4/63))*"PQ"))
Simulation type: assignment
+Funq(t)
Initial expression: (1-0.55)*("P680+/Qa-/Qb"+"P680+/Qa-/Qb-"+"P680+/Qa-/Qb2-"+"P680/Qa-/Qb2-"+"P680/Qa-/Qb"+"P680/Qa-/Qb-")/(1-0.55*("P680+/Qa-/Qb"+"P680+/Qa-/Qb-"+"P680+/Qa-/Qb2-"+"P680/Qa-/Qb2-"+"P680/Qa-/Qb"+"P680/Qa-/Qb-"))
Simulation type: assignment
+I820
Initial expression: 10^(-(2.16*10^(-7)*1590*"Pc+" + 2.16*10^(-7)*10300*("P700+/Fb"+"P700+/Fb-")))
Simulation type: assignment
