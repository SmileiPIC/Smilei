Sampling a Maxwell-Jüttner distribution
---------------------------------------

We base our method on that described in the Appendix B of
an `article by Schnittman and Krolik <http://dx.doi.org/10.1088/0004-637X/777/1/11>`_.

The Maxwell-Jüttner distribution, as a function of the Lorentz factor :math:`\gamma`, reads

.. math::
  
  f(\gamma) =  \gamma^2 \beta 
  \exp\left(- \frac {\gamma}{\theta} \right)

where :math:`\theta` is the temperature divided by :math:`mc^2`.
It is problematic that the change of variable :math:`\gamma/\theta` is impossible, because
it requires the cumulative distribution function to be computed for every different temperature.

Instead, the "rejection method" makes it possible to choose another function :math:`g(\gamma)`
such that :math:`g(\gamma)>f(\gamma)` everywhere. It can be chosen so that the cumulative
distribution function :math:`G(\gamma)` is easy to inverse. First, we take a random
number :math:`U_1` between 0 and 1, and sample the value :math:`\gamma_1=G^{-1}(U_1)`.
Second, we pick another random number :math:`U_2`, and if :math:`U_2<f(\gamma_1)/g(\gamma_1)`,
we keep the value :math:`\gamma_1`. Otherwise, we start over to choose another :math:`U_1`,
and so on until a good value is found.

In this particular case, we choose

.. math:: 
  
  g(\gamma) = \gamma^2 
  \exp\left(- \frac {\gamma}{\theta} \right)

which verifies :math:`g(\gamma)>f(\gamma)` and which has the cumulative distribution function

.. math::
  
  G(\gamma) = \int_1^\gamma g(x) dx = 1 - \exp\left[H(\gamma/\theta)-H(1/\theta)\right]
  
where :math:`H(u) = -u +\ln(1+u+u^2/2)`.

The rejection methods proceeds as

1. pick a random :math:`U_1`
2. calculate :math:`\gamma_1=G^{-1}(U_1)=\theta\; H^{-1}[\ln(1-U_1)+H(1/\theta)]`
3. pick a random :math:`U_2`
4. select :math:`\gamma_1` if :math:`U_2<\sqrt{1-\gamma_1^{-2}}`, otherwise restart from point 1

Now, to do this, we need to know :math:`H^{-1}`, which is not easy. We choose to tabulate it
in Smilei. For :math:`X>-\exp(-26)`, we use the series development :math:`H^{-1}(X) = (-6X)^{1/3}`.
For :math:`X<-\exp(12)`, we use the fit :math:`H^{-1}(X) = -X + 11.35(-X)^{0.06}`.
For all points in between, the function is linearly interpolated in log-log scale over 1000
tabulated values.

Note that the rejection method requires to pick several random numbers if the functions
:math:`f` and :math:`g` differ significantly. This strongly slows the calculation down
when the temperature is non-relativistic. For this reason, we fall back to the 
Maxwell-Boltzmann distribution when :math:`\theta<0.1`.