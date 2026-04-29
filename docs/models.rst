Fit models
==========

This page summarizes the analytical expressions implemented in
``TGAnalysis.model1``, ``TGAnalysis.model2``, and ``TGAnalysis.model3``.

Notation
--------

* :math:`t`: time (ps).
* :math:`t_0`: time-zero shift.
* :math:`\sigma`: Gaussian instrument-response width.
* :math:`k, k_1, k_2`: decay rates in ps :math:`^{-1}`.
* :math:`\tau = 1/k`: decay time (reported in fs in fit summaries).
* :math:`\mathrm{erf}(\cdot)`: error function.

Model 1 (``model1``)
--------------------

A mono-exponential decay (rate form) convolved with a Gaussian response plus
a step-like term:

.. math::

   f_1(t) =
   A_1 \exp\!\left(-(t-t_0)k\right)
   \exp\!\left(\frac{(k\sigma)^2}{2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k\sigma^2}{\sqrt{2}\sigma}\right)\right]
   +
   A_2 \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right].

Model 2 (``model2``)
--------------------

Model 1 plus a damped oscillatory contribution:

.. math::

   f_2(t) = f_{\mathrm{exp}}(t) + f_{\mathrm{step}}(t) + f_{\mathrm{osc}}(t),

where:

.. math::

   f_{\mathrm{exp}}(t) =
   A_1 \exp\!\left(-(t-t_0)k_1\right)
   \exp\!\left(\frac{(k_1\sigma)^2}{2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k_1\sigma^2}{\sqrt{2}\sigma}\right)\right],

.. math::

   f_{\mathrm{step}}(t) =
   A_2 \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right],

.. math::

   f_{\mathrm{osc}}(t) =
   A_3 \exp\!\left[-k_2(t-t_0)\right]
   \exp\!\left(-\frac{(k_2\sigma)^2}{2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k_2\sigma^2}{\sqrt{2}\sigma}\right)\right]
   \left[\sin(\omega t-\phi)-\cos(\omega t-\phi)\right].

Model 3 (``model3``)
--------------------

Extended multi-component form with two decay rates, two additional damped
oscillatory rates, and a Gaussian-convolved step term:

.. math::

   \text{with independent rates } k_1,\;k_2,\;k_{10},\;k_{20}

.. math::

   f_3(t) =
   A_1 e^{-(t-t_0)k_1} e^{-\frac{(k_1\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_1\sigma^2}{\sqrt{2}\sigma}\right)\right]
   +
   A_2 e^{-(t-t_0)k_2} e^{-\frac{(k_2\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_2\sigma^2}{\sqrt{2}\sigma}\right)\right]
   -
   A_3 e^{-(t-t_0)k_{10}} e^{-\frac{(k_{10}\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_{10}\sigma^2}{\sqrt{2}\sigma}\right)\right]\Psi(t)
   +
   A_3 e^{-(t-t_0)k_{20}} e^{-\frac{(k_{20}\sigma)^2}{2}}
   \left[1+\mathrm{erf}\!\left(\frac{t-t_0-k_{20}\sigma^2}{\sqrt{2}\sigma}\right)\right]\Psi(t)
   +
   A_4\left[1+\mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right],

.. math::

   \Psi(t)=\sin(\omega t-\phi)-\cos(\omega t-\phi).

Practical note
--------------

Use :meth:`tg_analysis.TGAnalysis.get_fit_parameters` to fit all loaded scans.
The method stores values and propagated errors in ``params_fit`` and reports
reduced :math:`\chi^2` and :math:`R^2` for each scan.
