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

   model_idxs = 2
   analysis.get_fit_parameters(
       model_idxs=model_idxs,
       initial_guess_bool=False,
       bounds=True,
   )

7. Plot fit results:

.. code-block:: python

   analysis.plot_fits()
   analysis.plot_taus_vs_energy()

8. Use stacked plots for detailed inspection:

.. code-block:: python

   analysis.plot_stacked_signals(
       limits_time=(0.0, 0.8),
       limits_signal=(0.1, 1.0),
       ylog_scale=False,
       plot_ind=[6],
       data_over_fit=False,
   )

What each block does in practice
--------------------------------

* ``TGAnalysis(json_path)`` loads metadata and builds ``df_scans`` with all ``ScanXXX`` entries.
* ``get_data_scan(...)`` reads ``.mat`` files from ``main_path``, applies optional time filtering, normalizes traces, and aligns the time origin.
* ``get_fit_parameters(...)`` runs ``curve_fit`` scan by scan and computes ``tau``, fit errors, reduced :math:`\chi^2`, and :math:`R^2`.
* ``plot_fits()`` compares measured data vs fitted curves and includes relative error.
* ``plot_taus_vs_energy()`` combines decay times with absorbance reference curves.
* ``plot_stacked_signals(...)`` compares experimental traces and fitted curves in stacked panels.

Recommended good practices
--------------------------

* Confirm that ``main_path`` in the JSON points to the correct ``.mat`` directory.
* Run ``get_data_scan`` before any fitting method or fit-related plot.
* Start with ``model_idxs = 2`` (as in the notebook), then compare with mixed model configurations.
* Keep exploratory work in notebooks and migrate stable steps to scripts/modules.
