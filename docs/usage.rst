Usage guide (based on results_tg.ipynb)
=======================================

This section documents the practical workflow used in ``src/results_tg.ipynb``.

Minimal workflow
----------------

1. Import the class:

.. code-block:: python

   from tg_analysis import TGAnalysis

2. Define analysis mode and target fitted parameter:

.. code-block:: python

   mode = "constant_I"  # or "constant_E"
   param_name = "tau"

3. Define the JSON metadata file:

.. code-block:: python

   json_path = "external_files/parameters_CoOx12.json"

4. Create the analysis object:

.. code-block:: python

   analysis = TGAnalysis(json_path)

5. Inspect phase-space coverage (energy vs intensity):

.. code-block:: python

   analysis.plot_phase_space(errors_bool=True)

6. Select scans according to mode:

.. code-block:: python

   if mode == "constant_E":
       params_scan = {"E": 63.6, "I": "all"}
   else:
       params_scan = {"E": "all", "I": 2.0}
   analysis.get_data_scan(params_scan)

7. Fit models and extract parameters:

.. code-block:: python

   model_idxs = 2
   analysis.get_fit_parameters(
       model_idxs=model_idxs,
       initial_guess_bool=True,
       bounds=True,
   )

8. Plot fit results:

.. code-block:: python

   analysis.plot_fits()
   if mode == "constant_E":
       analysis.plot_params_vs_intensity(param_name=param_name, errors_bool=False)
   else:
       analysis.plot_params_vs_energy(param_name=param_name, errors_bool=False)

9. Use stacked plots for detailed inspection:

.. code-block:: python

   analysis.plot_stacked_signals(
       limits_time=(0.0, 0.7),
       limits_signal=(0.01, 0.09),
       ylog_scale=True,
       plot_ind=[3],
       data_over_fit=True,
   )

10. Compare fit behavior across model configurations:

.. code-block:: python

   models_config = [
       {"model_idxs": 1, "initial_guess_bool": True, "bounds": True, "label_model": "Model 1", "color": "blue"},
       {"model_idxs": 2, "initial_guess_bool": True, "bounds": True, "label_model": "Model 2", "color": "red"},
       {"model_idxs": 3, "initial_guess_bool": True, "bounds": True, "label_model": "Model 3", "color": "green"},
   ]
   analysis.plot_params_all_models(
       models_config,
       param_name=param_name,
       mode=mode,
       errors_bool=False,
       save_path="../report/figures/example_all_models.png",
   )

``mode`` must match how scans were selected: ``"constant_I"`` (parameter vs energy) or ``"constant_E"`` (parameter vs intensity).

Batch export (notebook)
-------------------------

``results_tg.ipynb`` includes a loop over modes, models, and parameter names so you can regenerate report figures (fits, stacked signals, per-parameter plots, and optional ``plot_params_all_models`` overlays). Parameter names depend on the model (for example Model 1: ``amp1``, ``t0``, ``tau``, ``sigma``, ``ampoff``, ``r2``; Model 2: ``amp1``, ``t0``, ``tau``, ``sigma``, ``ampoff``, ``amp2``, ``tau2``, ``r2``; Model 3: ``amp1``, ``amp2``, ``amp3``, ``t0``, ``tau``, ``tau2``, ``tau3``, ``sigma``, ``r2``). Adjust paths and ``models_config`` filters if a parameter exists only for some models.

What each block does in practice
--------------------------------

* ``TGAnalysis(json_path)`` loads metadata and builds ``df_scans`` with all ``ScanXXX`` entries.
* ``get_data_scan(...)`` reads ``.mat`` files from ``main_path``, applies optional time filtering, subtracts baseline, and normalizes each TG trace by its integral.
* ``get_fit_parameters(...)`` runs ``curve_fit`` scan by scan and stores values + uncertainties in ``params_fit`` for each parameter.
* ``plot_fits()`` compares measured data vs fitted curves and reports :math:`\chi^2/\mathrm{dof}` and :math:`R^2`.
* ``plot_params_vs_energy(...)`` overlays absorbance references and any selected fitted parameter when using constant-intensity mode.
* ``plot_params_vs_intensity(...)`` plots fitted parameters against intensity when using constant-energy mode.
* ``plot_params_all_models(...)`` compares one fitted parameter across multiple model configurations with explicit scan mode.
* ``plot_stacked_signals(...)`` compares experimental traces and fitted curves in stacked panels with fit-centered time shifts.

Recommended good practices
--------------------------

* Confirm that ``main_path`` in the JSON points to the correct ``.mat`` directory.
* Run ``get_data_scan`` before any fitting method or fit-related plot.
* Start with one global model index (for example ``model_idxs = 2``), then compare alternatives using ``plot_params_all_models(..., mode=mode)``.
* Keep exploratory work in notebooks and migrate stable steps to scripts/modules.
