Usage guide (based on results_tg.ipynb)
=======================================

This section documents the practical workflow used in ``src/results_tg.ipynb``.

Minimal workflow
----------------

1. Import the class:

.. code-block:: python

   from tg_analysis import TGAnalysis

2. Define the JSON metadata file:

.. code-block:: python

   json_path = "external_files/parameters_CoOx12.json"

3. Create the analysis object:

.. code-block:: python

   analysis = TGAnalysis(json_path)

4. Inspect phase-space coverage (energy vs intensity):

.. code-block:: python

   analysis.plot_phase_space()

5. Select scans at fixed intensity (or fixed energy):

.. code-block:: python

   params_scan = {"E": "all", "I": 2.0}
   analysis.get_data_scan(params_scan)

6. Fit models and extract parameters:

.. code-block:: python

   model_idxs = 1
   analysis.get_fit_parameters(
       model_idxs=model_idxs,
       initial_guess_bool=True,
       bounds=True,
   )

7. Plot fit results:

.. code-block:: python

   analysis.plot_fits()
   analysis.plot_params_vs_energy(param_name="tau", errors_bool=True)

8. Use stacked plots for detailed inspection:

.. code-block:: python

   analysis.plot_stacked_signals(
       limits_time=(0.0, 0.8),
       limits_signal=(0.1, 1.0),
       ylog_scale=False,
       plot_ind=None,
       data_over_fit=True,
   )

9. Compare fit behavior across model configurations:

.. code-block:: python

   models_config = [
       {"model_idxs": 1, "initial_guess_bool": True, "bounds": True, "label_model": "Model 1"},
       {"model_idxs": 2, "initial_guess_bool": True, "bounds": True, "label_model": "Model 2"},
   ]
   analysis.plot_params_all_models(models_config, param_name="tau", errors_bool=False)

What each block does in practice
--------------------------------

* ``TGAnalysis(json_path)`` loads metadata and builds ``df_scans`` with all ``ScanXXX`` entries.
* ``get_data_scan(...)`` reads ``.mat`` files from ``main_path``, applies optional time filtering, subtracts baseline, and normalizes each TG trace by its integral.
* ``get_fit_parameters(...)`` runs ``curve_fit`` scan by scan and stores values + uncertainties in ``params_fit`` for each parameter.
* ``plot_fits()`` compares measured data vs fitted curves and reports :math:`\chi^2/\mathrm{dof}` and :math:`R^2`.
* ``plot_params_vs_energy(...)`` overlays absorbance references and any selected fitted parameter (for example ``"tau"``, ``"sigma"``, or ``"omega"``).
* ``plot_params_all_models(...)`` compares one fitted parameter across multiple model configurations.
* ``plot_stacked_signals(...)`` compares experimental traces and fitted curves in stacked panels.

Recommended good practices
--------------------------

* Confirm that ``main_path`` in the JSON points to the correct ``.mat`` directory.
* Run ``get_data_scan`` before any fitting method or fit-related plot.
* Start with one global model index (for example ``model_idxs = 1``), then compare alternatives using ``plot_params_all_models(...)``.
* Keep exploratory work in notebooks and migrate stable steps to scripts/modules.
