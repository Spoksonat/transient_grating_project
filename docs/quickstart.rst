Quickstart
==========

Use this page to run a first end-to-end analysis in a few minutes.

Prerequisites
-------------

* The project dependencies are installed.
* ``external_files/parameters_CoOx12.json`` exists and its ``main_path`` points to valid ``.mat`` files.

Minimal example
---------------

.. code-block:: python

   from tg_analysis import TGAnalysis

   json_path = "external_files/parameters_CoOx12.json"
   analysis = TGAnalysis(json_path)

   # 1) Explore available scan conditions
   analysis.plot_phase_space()

   # 2) Load traces at fixed intensity
   params_scan = {"E": "all", "I": 2.0}
   analysis.get_data_scan(params_scan)

   # 3) Fit model 2 and inspect results
   analysis.get_fit_parameters(model_idxs=2, initial_guess_bool=False, bounds=True)
   analysis.plot_fits()
   analysis.plot_taus_vs_energy()

Common next steps
-----------------

* Try ``params_scan = {"E": value, "I": "all"}`` to compare intensities at fixed energy.
* Compare models with ``plot_taus_all_models(...)``.
* Use ``plot_stacked_signals(...)`` to inspect specific scans.
