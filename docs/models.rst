Fit models
==========

This page summarizes the analytical expressions implemented in
``TGAnalysis.model1``, ``TGAnalysis.model2``, and ``TGAnalysis.model3``.

Notation
--------

* :math:`t`: time (ps).
* :math:`t_0`: time-zero shift.
* :math:`\sigma`: Gaussian instrument-response width.
* :math:`k, k_1, k_2, k_3`: decay rates in ps :math:`^{-1}`.
* :math:`\tau = 1/k`: decay time (reported in fs in fit summaries).
* :math:`\mathrm{erf}(\cdot)`: error function.

Model 1 (``model1``)
--------------------

A mono-exponential decay convolved with a Gaussian response plus
an offset step term:

.. math::

   f_1(t) =
   A_1 \exp\!\left(-(t-t_0)k\right)
   \exp\!\left(\frac{(k\sigma)^2}{2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k\sigma^2}{\sqrt{2}\sigma}\right)\right]
   +
   A_{\mathrm{off}} \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right].

Model 2 (``model2``)
--------------------

Model 1 plus a second exponential decay channel:

.. math::

   f_2(t) = f_{\mathrm{exp},1}(t) + f_{\mathrm{off}}(t) + f_{\mathrm{exp},2}(t),

where:

.. math::

   f_{\mathrm{exp},1}(t) =
   A_1 \exp\!\left(-(t-t_0)k_1\right)
   \exp\!\left(\frac{(k_1\sigma)^2}{2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k_1\sigma^2}{\sqrt{2}\sigma}\right)\right],

.. math::

   f_{\mathrm{off}}(t) =
   A_{\mathrm{off}} \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right],

.. math::

   f_{\mathrm{exp},2}(t) =
   A_2 \exp\!\left[-k_2(t-t_0)\right]
   \exp\!\left(-\frac{(k_2\sigma)^2}{2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k_2\sigma^2}{\sqrt{2}\sigma}\right)\right].

Model 3 (``model3``)
--------------------

A sum of three Gaussian-convolved exponential channels:

.. math::

   f_3(t) =
   A_1 e^{-(t-t_0)k_1} e^{-\frac{(k_1\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_1\sigma^2}{\sqrt{2}\sigma}\right)\right]
   +
   A_2 e^{-(t-t_0)k_2} e^{-\frac{(k_2\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_2\sigma^2}{\sqrt{2}\sigma}\right)\right]
   +
   A_3 e^{-(t-t_0)k_3} e^{-\frac{(k_3\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_3\sigma^2}{\sqrt{2}\sigma}\right)\right].

Bounds and reported :math:`\tau`
---------------------------------

Fitted decay rates :math:`k` (in ps :math:`^{-1}`) are converted to
:math:`\tau = 1000/k` fs when filling ``"tau"``, ``"tau2"``, and ``"tau3"``
entries in ``params_fit``. Model 1 and Model 2 also expose ``"ampoff"``
for the convolved offset term; for plotting compatibility across models,
Model 3 stores ``"ampoff"`` as ``NaN``.

Practical note
--------------

Use :meth:`tg_analysis.TGAnalysis.get_fit_parameters` to fit all loaded scans.
The method stores values and propagated errors in ``params_fit`` and reports
reduced :math:`\chi^2` and :math:`R^2` for each scan.
