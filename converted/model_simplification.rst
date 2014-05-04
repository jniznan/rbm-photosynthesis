
.. code:: python

    from models.holzwarth.naive import model
    from tools import dfs, rules
    from pysb import bng

.. parsed-literal::

    /usr/lib64/python2.7/site-packages/mpl_toolkits/__init__.py:2: UserWarning: Module argparse was already imported from /usr/lib64/python2.7/argparse.pyc, but /usr/lib/python2.7/site-packages is being added to sys.path
      __import__('pkg_resources').declare_namespace(__name__)


Original model
--------------

.. code:: python

    rules(model)



.. parsed-literal::

    [Rule('deexcitation_ChlD1_fluorescence', PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='exc') >> PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='n'), deex_flu),
     Rule('deexcitation_ChlD1_heat', PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='exc') >> PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='n'), deex_heat),
     Rule('deexcitation_antenna_fluorescence', PSII(P680='n', Qa='n', PhD1='n', ac='exc', ChlD1='n') >> PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='n'), deex_flu),
     Rule('deexcitation_antenna_heat', PSII(P680='n', Qa='n', PhD1='n', ac='exc', ChlD1='n') >> PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='n'), deex_heat),
     Rule('excitation_transfer', PSII(P680='n', Qa='n', PhD1='n', ac='exc', ChlD1='n') <> PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='exc'), extrans, extrans_re),
     Rule('primary_charge_separation_recombination', PSII(P680='n', Qa='n', PhD1='n', ac='n', ChlD1='exc') <> PSII(P680='n', Qa='n', PhD1='m', ac='n', ChlD1='p'), ch_se, ch_re),
     Rule('quinone_qa_reduction_oxidation', PSII(P680='p', Qa='n', PhD1='m', ac='n', ChlD1='n') <> PSII(P680='p', Qa='m', PhD1='n', ac='n', ChlD1='n'), qa_red, qa_oxi),
     Rule('stable_pair_generation_degeneration', PSII(P680='n', Qa='n', PhD1='m', ac='n', ChlD1='p') <> PSII(P680='p', Qa='n', PhD1='m', ac='n', ChlD1='n'), stab_gen, stab_deg)]



Simplified model
----------------

.. code:: python

    m1 = dfs(model)
.. code:: python

    rules(m1)



.. parsed-literal::

    [Rule('deexcitation_ChlD1_fluorescence', PSII(ChlD1='exc') >> PSII(ChlD1='n'), deex_flu),
     Rule('deexcitation_ChlD1_heat', PSII(ChlD1='exc') >> PSII(ChlD1='n'), deex_heat),
     Rule('deexcitation_antenna_fluorescence', PSII(ac='exc') >> PSII(ac='n'), deex_flu),
     Rule('deexcitation_antenna_heat', PSII(ac='exc') >> PSII(ac='n'), deex_heat),
     Rule('excitation_transfer', PSII(ac='exc', ChlD1='n') <> PSII(ac='n', ChlD1='exc'), extrans, extrans_re),
     Rule('primary_charge_separation_recombination', PSII(PhD1='n', ChlD1='exc') <> PSII(PhD1='m', ChlD1='p'), ch_se, ch_re),
     Rule('quinone_qa_reduction_oxidation', PSII(Qa='n', PhD1='m', ChlD1='n') <> PSII(Qa='m', PhD1='n', ChlD1='n'), qa_red, qa_oxi),
     Rule('stable_pair_generation_degeneration', PSII(P680='n', PhD1='m', ChlD1='p') <> PSII(P680='p', PhD1='m', ChlD1='n'), stab_gen, stab_deg)]



BNGL original
-------------

.. code:: python

    print bng.BngGenerator(model).get_content()

.. parsed-literal::

    begin parameters
      deex_heat    2.000000e+05
      qa_oxi       2.400000e+09
      ch_se        5.000000e+11
      stab_gen     2.500000e+11
      qa_red       4.800000e+09
      extrans      1.920000e+10
      extrans_re   2.500000e+10
      deex_flu     6.700000e+07
      ch_re        5.000000e+11
      stab_deg     6.700000e+10
      init         1.000000e+00
    end parameters
    
    begin molecule types
      PSII(P680~n~p,Qa~n~m,PhD1~n~m,ac~n~exc,ChlD1~n~p~exc)
    end molecule types
    
    begin observables
      Molecules Ant_exc     PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n)
      Molecules ChlD1_exc   PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~exc)
    end observables
    
    begin species
      PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n)   init
    end species
    begin reaction rules
      deexcitation_antenna_fluorescence:        PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n) -> PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~n)    deex_flu
      deexcitation_antenna_heat:                PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n) -> PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~n)    deex_heat
      deexcitation_ChlD1_fluorescence:          PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~exc) -> PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~n)    deex_flu
      deexcitation_ChlD1_heat:                  PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~exc) -> PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~n)    deex_heat
      excitation_transfer:                      PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n) <-> PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~exc)    extrans, extrans_re
      primary_charge_separation_recombination:  PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~exc) <-> PSII(P680~n,Qa~n,PhD1~m,ac~n,ChlD1~p)    ch_se, ch_re
      stable_pair_generation_degeneration:      PSII(P680~n,Qa~n,PhD1~m,ac~n,ChlD1~p) <-> PSII(P680~p,Qa~n,PhD1~m,ac~n,ChlD1~n)    stab_gen, stab_deg
      quinone_qa_reduction_oxidation:           PSII(P680~p,Qa~n,PhD1~m,ac~n,ChlD1~n) <-> PSII(P680~p,Qa~m,PhD1~n,ac~n,ChlD1~n)    qa_red, qa_oxi
    end reaction rules
    
    


BNGL simplified
---------------

.. code:: python

    print bng.BngGenerator(m1).get_content()

.. parsed-literal::

    begin parameters
      deex_heat    2.000000e+05
      qa_oxi       2.400000e+09
      ch_se        5.000000e+11
      stab_gen     2.500000e+11
      qa_red       4.800000e+09
      extrans      1.920000e+10
      extrans_re   2.500000e+10
      deex_flu     6.700000e+07
      ch_re        5.000000e+11
      stab_deg     6.700000e+10
      init         1.000000e+00
    end parameters
    
    begin molecule types
      PSII(P680~n~p,Qa~n~m,PhD1~n~m,ac~n~exc,ChlD1~n~p~exc)
    end molecule types
    
    begin observables
      Molecules Ant_exc     PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n)
      Molecules ChlD1_exc   PSII(P680~n,Qa~n,PhD1~n,ac~n,ChlD1~exc)
    end observables
    
    begin species
      PSII(P680~n,Qa~n,PhD1~n,ac~exc,ChlD1~n)   init
    end species
    begin reaction rules
      stable_pair_generation_degeneration:      PSII(P680~n,PhD1~m,ChlD1~p) <-> PSII(P680~p,PhD1~m,ChlD1~n)    stab_gen, stab_deg
      quinone_qa_reduction_oxidation:           PSII(Qa~n,PhD1~m,ChlD1~n) <-> PSII(Qa~m,PhD1~n,ChlD1~n)    qa_red, qa_oxi
      excitation_transfer:                      PSII(ac~exc,ChlD1~n) <-> PSII(ac~n,ChlD1~exc)    extrans, extrans_re
      primary_charge_separation_recombination:  PSII(PhD1~n,ChlD1~exc) <-> PSII(PhD1~m,ChlD1~p)    ch_se, ch_re
      deexcitation_antenna_fluorescence:        PSII(ac~exc) -> PSII(ac~n)    deex_flu
      deexcitation_antenna_heat:                PSII(ac~exc) -> PSII(ac~n)    deex_heat
      deexcitation_ChlD1_fluorescence:          PSII(ChlD1~exc) -> PSII(ChlD1~n)    deex_flu
      deexcitation_ChlD1_heat:                  PSII(ChlD1~exc) -> PSII(ChlD1~n)    deex_heat
    end reaction rules
    
    


