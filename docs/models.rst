Fit models
==========

This page summarizes the analytical expressions implemented in
``TGAnalysis.model1``, ``TGAnalysis.model2``, and ``TGAnalysis.model3``.

Notation
--------

* :math:`t`: time (ps).
* :math:`t_0`: time-zero shift.
* :math:`\sigma`: Gaussian instrument-response width.
* :math:`\tau, \tau_2`: decay constants.
* :math:`\mathrm{erf}(\cdot)`: error function.

Model 1 (``model1``)
--------------------

A mono-exponential decay convolved with a Gaussian response plus a step-like
background term:

.. math::

   f_1(t) =
   A_1 \exp\!\left(-\frac{t-t_0}{\tau}\right)
   \exp\!\left(\frac{\sigma^2}{2\tau^2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-\sigma^2/\tau}{\sqrt{2}\sigma}\right)\right]
   +
   A_2 \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right].

Model 2 (``model2``)
--------------------

Model 1 plus a damped oscillatory term:

.. math::

   f_2(t) = f_{\mathrm{exp}}(t) + f_{\mathrm{step}}(t) + f_{\mathrm{osc}}(t),

where:

.. math::

   f_{\mathrm{exp}}(t) =
   A_1 \exp\!\left(-\frac{t-t_0}{\tau}\right)
   \exp\!\left(\frac{\sigma^2}{2\tau^2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-\sigma^2/\tau}{\sqrt{2}\sigma}\right)\right],

.. math::

   f_{\mathrm{step}}(t) =
   A_2 \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right],

.. math::

   f_{\mathrm{osc}}(t) =
   A_3 \exp\!\left[-k(t-t_0)\right]
   \exp\!\left(-\frac{1}{2}k\sigma^2\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-k\sigma^2}{\sqrt{2}\sigma}\right)\right]
   \left[\sin(\omega t-\phi)-\cos(\omega t-\phi)\right].

Model 3 (``model3``)
--------------------

A bi-exponential convolved form with a shared Gaussian response:

.. math::

   f_3(t) =
   A_1 \left[1 + \mathrm{erf}\!\left(\frac{t-t_0}{\sqrt{2}\sigma}\right)\right]
   +
   A_2 \exp\!\left(-\frac{t-t_0}{\tau}\right)
   \exp\!\left(\frac{\sigma^2}{2\tau^2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-\sigma^2/\tau}{\sqrt{2}\sigma}\right)\right]
   +
   A_3 \exp\!\left(-\frac{t-t_0}{\tau_2}\right)
   \exp\!\left(\frac{\sigma^2}{2\tau_2^2}\right)
   \left[1 + \mathrm{erf}\!\left(\frac{t-t_0-\sigma^2/\tau_2}{\sqrt{2}\sigma}\right)\right].

Practical note
--------------

Use :meth:`tg_analysis.TGAnalysis.get_fit_parameters` to apply these models to
all loaded scans and compute fitted parameters, uncertainties, reduced
:math:`\chi^2` values, and :math:`R^2` values.
