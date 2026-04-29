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

   mode = "constant_I"
   param_name = "tau"

   json_path = "external_files/parameters_CoOx12.json"
   analysis = TGAnalysis(json_path)

   # 1) Explore available scan conditions
   analysis.plot_phase_space(errors_bool=True)

   # 2) Load traces according to selected mode
   if mode == "constant_E":
       params_scan = {"E": 63.6, "I": "all"}
   else:
       params_scan = {"E": "all", "I": 2.0}
   analysis.get_data_scan(params_scan)

   # 3) Fit model and inspect results
   analysis.get_fit_parameters(model_idxs=2, initial_guess_bool=True, bounds=True)
   analysis.plot_fits()
   if mode == "constant_E":
       analysis.plot_params_vs_intensity(param_name=param_name, errors_bool=False)
   else:
       analysis.plot_params_vs_energy(param_name=param_name, errors_bool=False)

Common next steps
-----------------

* Try ``params_scan = {"E": value, "I": "all"}`` to compare intensities at fixed energy.
* Compare models with ``plot_params_all_models(..., mode=mode)``.
* Use ``plot_stacked_signals(...)`` to inspect specific scans.
