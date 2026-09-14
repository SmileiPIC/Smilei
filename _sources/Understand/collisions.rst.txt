Binary collisions & reactions
-----------------------------


Relativistic binary collisions between particles have been implemented in
:program:`Smilei` following these references:

* [Perez2012]_: overview of the technique.
* [Nanbu1997]_ and [Nanbu1998]_: the original approach.
* [Sentoku2008]_, [Lee1984]_ and [Frankel1979]_: additional information.
* The correction suggested in [Higginson2020]_ has been applied since v4.5.
* [to be published] A correction for the screening of bound electrons.
  This makes a significant difference for weakly-ionized atoms.

This collision scheme can host reactions between the colliding
macro-particles, when requested:

* Ionization of an atom by collision with an electron.
* Nuclear reaction between two atoms.

Please refer to :ref:`that doc <Collisions>` for an explanation of how to add
collisions in the namelist file.


----

The binary collision scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Collisions are calculated at each timestep and for each collision block
given in the input file:

* Macro-particles that should collide are randomly paired.
* Average parameters are calculated (densities, etc.).
* For each pair:
  
  * Calculate the momenta in the center-of-mass (COM) frame.
  * Calculate the coulomb log if requested (see [Perez2012]_).
  * If the collision corresponds to a nuclear reaction (optional),
    the reaction probability is computed and new particles are created
    if successful.
  * Calculate the collision rate.
  * Randomly pick the deflection angle.
  * Deflect particles in the COM frame and switch back to the laboratory frame.
  * If the collision corresponds to ionization (optional),
    its probability is computed and new electrons are created
    if successful.

.. rubric:: Modifications in Smilei

* A typo from [Perez2012]_ is corrected: in Eq. (22), corresponding to
  the calculation of the Coulomb Logarithm, the last parenthesis is
  written as a squared expression, but should not.

* The deflection angle distribution given by [Nanbu1997]_
  (which is basically a fit from Monte-Carlo simulations)
  is modified for better accuracy and performance.
  Given Nanbu's :math:`s` parameter and a random number :math:`U\in [0,1]`,
  the deflection angle :math:`\chi` is:
  
  .. math::
  
    \sin^2\frac\chi 2 = \begin{cases} 
    \alpha U/\sqrt{1-U + \alpha^2 U} &, s < 4\\
    1-U &, \textrm{otherwise}
    \end{cases}
    
  where :math:`\alpha = 0.37 s-0.005 s^2-0.0064 s^3`.

----

Test cases for collisions
^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: 1. Beam relaxation

An electron beam with narrow energy spread enters an ion background with :math:`T_i=10` eV.
The ions are of very small mass :math:`m_i=10 m_e` to speed-up the calculation.
Only e-i collisions are calculated.
The beam gets strong isotropization => the average velocity relaxes to zero.

Three figures show the time-evolution of the longitudinal :math:`\left<v_\|\right>`
and transverse velocity :math:`\sqrt{\left<v_\perp^2\right>}`

* :numref:`beam1` : initial velocity = 0.05, ion charge = 1
* :numref:`beam2` : initial velocity = 0.01, ion charge = 1
* :numref:`beam3` : initial velocity = 0.01, ion charge = 3

Each of these figures show 3 different blue and red curves which correspond to different
ratios of particle weights: 0.1, 1, and 10.

.. _beam1:

.. figure:: /_static/beam_relaxation123.png
  :width: 10cm
  
  Relaxation of an electron beam. Initial velocity = 0.05, ion charge = 1.
  
.. _beam2:

.. figure:: /_static/beam_relaxation456.png
  :width: 10cm
  
  Relaxation of an electron beam. Initial velocity = 0.01, ion charge = 1.

.. _beam3:

.. figure:: /_static/beam_relaxation789.png
  :width: 10cm
  
  Relaxation of an electron beam. Initial velocity = 0.01, ion charge = 3.


The black lines correspond to the theoretical rates taken from the NRL formulary:

.. math::
  
  \nu_\| = -\left(1+\frac{m_e}{m_i}\right)\nu_0
  \quad\textrm{and}\quad
  \nu_\perp = 2\;\nu_0
  \quad\textrm{where}\quad
  \nu_0=\frac{e^4\,Z^{\star 2}\,n_i\,\ln\Lambda } { 4 \pi \epsilon_0^2 \,m_e^2\,v_e^3 }


The distribution is quickly non-Maxwellian so that theory is valid only at the beginning.


.. rubric:: 2. Thermalization

A population of electrons has a different temperature from that of the ion population.
Through e-i collisions, the two temperatures become equal.
The ions are of very small mass :math:`m_i=10 m_e` to speed-up the calculation.
Three cases are simulated, corresponding to different ratios of weights: 0.2, 1 and 5.
They are plotted in :numref:`thermalization`.

.. _thermalization:

.. figure:: /_static/thermalisation_ei123.png
  :width: 9cm
  
  Thermalization between two species.

The black lines correspond to the theoretical rates taken from the NRL formulary:

.. math::
  
  \nu_\epsilon = \frac{2}{3}\sqrt\frac{2}{\pi}
  \frac{e^4\,Z^{\star 2} \sqrt{m_em_i}\,n_i\,\ln\Lambda }
  { 4 \pi\epsilon_0^2 \,\left(m_eT_e+m_iT_i\right)^{3/2} }




.. rubric:: 3. Temperature isotropization

Electrons have a longitudinal temperature different from their transverse temperature.
They collide only with themselves (intra-collisions) and the anisotropy disappears
as shown in :numref:`temperature_isotropization`.

.. _temperature_isotropization:

.. figure:: /_static/temperature_isotropization1.png
  :width: 10cm
  
  Temperature isotropization of an electron population.

The black lines correspond to the theoretical rates taken from the NRL formulary:

.. math::
  
  \nu_T=\frac{e^4 \,n_e\,\ln\Lambda } { 8\pi^{3/2} \epsilon_0^2 \,m_e^{1/2}T_\|^{3/2} }
  A^{-2} \left[-3+(3-A)\frac{\rm{arctanh}(\sqrt{A})}{\sqrt{A}}\right]
  \quad \rm{where}\quad A=1-\frac{T_\perp}{T_\|}



.. rubric:: 4. Maxwellianization

Electrons start with zero temperature along :math:`y` and :math:`z`.
Their velocity distribution along :math:`x` is rectangular.
They collide only with themselves and the rectangle becomes a maxwellian 
as shown in :numref:`maxwellianization`.

.. _maxwellianization:

.. figure:: /_static/Maxwellianization1.png
  :width: 10cm
  
  Maxwellianization of an electron population.
  Each blue curve is the distribution at a given time.
  The red curve is an example of a gaussian function.



.. rubric:: 5. Stopping power

Test electrons (very low density) collide with background electrons of density
:math:`10\,n_c` and :math:`T_e=5` keV.
Depending on their initial velocity, they are slowed down at different rates,
as shown in :numref:`stoppingpower`.

.. _stoppingpower:

.. figure:: /_static/Stopping_power123.png
  :width: 10cm
  
  Stopping power of test electrons into a background electron population.
  Each point is one simulation. The black line is Frankel's theory [Frankel1979]_.


.. rubric:: 6. Conductivity

Solid-density Cu is simulated at different temperatures (e-i equilibrium) with only
e-i collisions. An electric field of :math:`E=3.2` GV/m (0.001 in code units) is
applied using two charged layers on each side of the solid Cu.
The electron velocity increases until a limit value :math:`v_f`.
The resulting conductivity :math:`\sigma=en_ev_f/E` is compared in
:numref:`conductivity` to the models in [Lee1984]_ and [Perez2012]_.

.. _conductivity:

.. figure:: /_static/conductivity.png
  :width: 10cm
  
  Conductivity of solid-density copper. Each point is one simulation.


----

.. _CollIonization:

Collisional ionization
^^^^^^^^^^^^^^^^^^^^^^

The binary collisions can also be ionizing if they are **electron-ion** collisions.
The approach is almost the same as that provided in [Perez2012]_.

When ionization is requested by setting ``ionizing=True``, a few additional operations
are executed:

* At the beginning of the run, cross-sections are calculated from tabulated binding
  energies (available for ions up to atomic number 100). These cross-sections are then
  tabulated for each requested ion species.
* Each timestep, the particle density :math:`n = n_e n_i/n_{ei}`
  (similar to the densities above for collisions) is calculated.
* During each collision, a probability for ionization is computed. If successful, 
  the ion charge is increased, the incident electron is slowed down, and a new electron
  is created.

.. rubric:: Warnings

* This scheme does not account for recombination, which would balance ionization
  over long time scales.

.. rubric:: Relativistic change of frame

A modification has been added to the theory of [Perez2012]_ in order to account for the
laboratory frame being different from the ion frame. Considering :math:`\overrightarrow{p_e}`
and :math:`\overrightarrow{p_i}` the electron and ion momenta in the laboratory frame, 
and their associated Lorentz factors :math:`\gamma_e` and :math:`\gamma_i`, we define
:math:`\overrightarrow{q_e}=\overrightarrow{p_e}/(m_e c)` and
:math:`\overrightarrow{q_i}=\overrightarrow{p_i}/(m_i c)`.
The Lorentz factor of the electron in the ion frame is 
:math:`\gamma_e^\star=\gamma_e\gamma_i-\overrightarrow{q_e}\cdot\overrightarrow{q_i}`.
The probability for ionization reads:

.. math::
  
  P = 1-\exp\left( - v_e \sigma n \Delta t \right) = 1-\exp\left( -V^\star \sigma^\star n \Delta t \right)

where :math:`v_e` is the electron velocity in the laboratory frame,
:math:`\sigma` is the cross-section in the laboratory frame, :math:`\sigma^\star`
is the cross-section in the ion frame, and 
:math:`V^\star=c\sqrt{\gamma_e^{\star\,2}-1}/(\gamma_e\gamma_i)`.

The loss of energy :math:`m_ec^2 \delta\gamma` of the incident electron translates into a change in momentum
:math:`{q_e^\star}' = \alpha_e q_e^\star` in the ion frame, with
:math:`\alpha_e=\sqrt{(\gamma_e^\star-\delta\gamma)^2-1}/\sqrt{\gamma_e^{\star2}-1}`.
In the laboratory frame, it becomes
:math:`\overrightarrow{q_e'}=\alpha_e\overrightarrow{q_e}+((1-\alpha_e)\gamma_e^\star-\delta\gamma)\overrightarrow{q_i}`.

A similar operation is done for defining the momentum of the new electron in the lab frame.
It is created with energy :math:`m_ec^2 (\gamma_w-1)` and its momentum is
:math:`q_w^\star = \alpha_w q_e^\star` in the ion frame, with
:math:`\alpha_w=\sqrt{\gamma_w^2-1}/\sqrt{\gamma_e^{\star 2}-1}`.
In the laboratory frame, it becomes
:math:`\overrightarrow{q_w}=\alpha_w\overrightarrow{q_e}+(\gamma_w-\alpha_w\gamma_e^\star)\overrightarrow{q_i}`.


.. rubric:: Multiple ionization

A modification has been added to the theory of [Perez2012]_ in order to account for 
multiple ionization in a single timestep. The approach for field ionization
by `Nuter et al <http://dx.doi.org/10.1063/1.3559494>`_
has been adapted to calculate the successive impact ionization probabilities
when an ion is ionized several times in a row.

Writing the probability to not ionize an ion already ionized :math:`i` times as
:math:`\bar{P}^i = \exp\left( -W_i\Delta t\right)`, and defining 
:math:`R^m_n = (1-W_m/W_n)^{-1}`, we can calculate the probability to ionize :math:`k` times
the ion:

.. math::
  
  P^i_k = \left\{
  \begin{array}{ll}
  \bar{P}^i
  &
  \quad\mathrm{if}\quad k=0
  \\
  \sum\limits_{p=0}^{k-1} R^{i+k}_{i+p} \left(\bar{P}^{i+k} - \bar{P}^{i+p}\right)
  \prod\limits_{j=0,j\ne p}^{k-1} R^{i+p}_{i+j}
  &
  \quad\mathrm{if}\quad 0<k<Z-Z^\star
  \end{array}
  \right.

..
  \\
  \sum\limits_{p=0}^{k-1} \left[ 1+R^{i+k}_{i+p}\left(\frac{W_{i+k}}{W_{i+p}}\bar{P}^{i+p} - \bar{P}^{i+k}\right) \right]
  \prod\limits_{j=0,j\ne p}^{k-1} R^{i+p}_{i+j}
  &
  \quad\mathrm{if}\quad k=k_\mathrm{max}


To simplify the calculation of :math:`P^i_k` (in particular the second case in the
equation above) we use the following equivalent expression:

.. math::
  
  P^i_k = A_{k-1} \sum\limits_{p=0}^{k-1}  \left(\bar{P}^{i+k} - \bar{P}^{i+p}\right)
  / B_k^p
  \quad\mathrm{if}\quad 0<k<Z-Z^\star

where :math:`A_k = \prod\limits_{j=0}^{k} W_{i+j}` and
:math:`B_k^p = \prod\limits_{j=0,j\ne p}^{k} (W_{i+j}-W_{i+p})`.
These two quantities can be computed recursively for each :math:`k`.


The cumulative probability :math:`F^i_k = \sum_{j=0}^{k} P^i_j` provides an efficient
way to pick when the ionization stops: we pick a random number :math:`U\in [0,1]` and
loop from :math:`k=0` to :math:`k_\mathrm{max}`. We stop ionizing when :math:`F^i_k>U`.



----

Test cases for ionization
^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: 1. Ionization rate

A cold plasma of :math:`\mathrm{Al}^{3+}` is set with density :math:`n_e=10^{21} \mathrm{cm}^{-3}`
and with all electrons drifting at a velocity :math:`v_e=0.03\,c`. The charge state of ions
versus time is shown in :numref:`IonizationRate` where the three dotted curves correspond
to three different weight ratios between electrons and ions.

.. _IonizationRate:

.. figure:: /_static/ionization_rate.png
  :width: 10cm
  
  Ionization of an aluminium plasma by drifting electrons.
  
The theoretical curve (in black) corresponds to :math:`1-\exp\left(v_en_e\sigma t\right)`
where :math:`\sigma` is the ionization cross section of :math:`\mathrm{Al}^{3+}` at the
right electron energy. The discrepancy at late time is due to the changing velocity
distributions and to the next level starting to ionize.


.. rubric:: 2. Inelastic stopping power

A cold, non-ionized Al plasma is set with density :math:`n_e=10^{21} \mathrm{cm}^{-3}`.
Electrons of various initial velocities are slowed down by ionizing collisions and their
energy loss is recorded as a function of time.

A few examples are given in the left graph of :numref:`IonizationStoppinPower`.
The theoretical curve is obtained from [Rohrlich1954]_. Note that this theory does not
work below a certain average ionization energy, in our case :math:`\sim 200` eV.

.. _IonizationStoppinPower:

.. figure:: /_static/ionization_stopping_power.png
  :width: 14cm
  
  Left: ionization slowing down versus time, for electrons injected at various
  initial energies into cold Al. Right: corresponding stopping power versus initial
  electron energy.
  
In the same figure, the graph on the right-hand-side provides the stopping power value
in the same context, at different electron energies. It is compared to the same theory.


.. rubric:: 3. Multiple ionization

If the timestep is large, multiple ionization can occur, especially with cold high-Z
material and high-energy electrons. The multiple ionization algorithm is not perfect,
as it does not shuffle the particles for each ionization. Thus, good statistical
sampling is reached after several timesteps. To test the potential error,
we ran simulations of electrons at 1 MeV incident on cold atoms. The evolution of the
secondary electron density is monitored versus time in :numref:`IonizationMultiple`.

.. _IonizationMultiple:

.. figure:: /_static/ionization_multiple.png
  :width: 10cm
  
  Secondary electron density *vs* time, for cold plasmas traversed by a 1 MeV electron beam.

The solid lines correspond to a very-well resolved ionization, whereas the dashed lines
correspond to a large timestep. A difference is visible initially, but decreases
quickly as the statistical sampling increases and as the subsequent ionization
cross-sections decrease.


.. rubric:: 3. Effect of neglecting recombination

As recombination is not accounted for, we can expect excess ionization to occur indefinitely
without being balanced to equilibrium. However, in many cases, the recombination rate
is small and can be neglected over the duration of the simulation. We provide an example
that is relevant to picosecond-scale laser-plasma interaction. Plasmas initially at
a density of 10 times the critical density are given various initial temperatures.
Ionization initially increases while the temperature decreases, until, after a while,
their charge state stagnates (it still increases, but very slowly).
In :numref:`IonizationRecombination`, these results are compared to a Thomas-Fermi model
from [Desjarlais2001]_.

.. _IonizationRecombination:

.. figure:: /_static/ionization_recombination.png
  :width: 12cm
  
  Final charge state of various plasmas at various temperatures.

The model does not account for detailed ionization potentials. It provides a rough
approximation, and is particularly questionable for low temperatures or high Z.
We observe that Smilei's approach for impact ionization provides decent estimates
of the ionization state. Detailed comparison to atomic codes has not been done yet.

----

.. _BoundElectronScreening:

Bound electron screening
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Coulomb potential that is assumed in the collision theory of Nanbu does not
correctly apply to neutral atoms or weakly-ionized ions. Indeed, in addition
to the screening caused by free electrons (Debye screening), the potential
must account for the bound-electron screening.

In Smilei, a simple model is available in the case of electron-ion collisions.
It relies on the Thomas-Fermi characteristic length
:math:`\lambda_{TF} = (9\pi^2/16Z)^{1/3}a_0`,
where :math:`a_0` is the Bohr radius, to describe the extent of the screening
aroud the nucleus. With a few approximations, the model suggests to replace,
in the collisional frequency, the term :math:`Z^{\star 2}\ln\Lambda_D` by

.. math::
  
  Z^2\ln\Lambda_{TFD} = Z^{\star 2} \ln\Lambda_D + (Z^2-Z^{\star 2})\ln\Lambda_{TF}

where :math:`\ln\Lambda_D` is the usual (Debye-screening) Coulomb logarithm,
and :math:`\Lambda_{TF} = \lambda_{TF}\Lambda_D /\lambda_D`. This approach is
basic but captures the essential trend. Little data for comparison is available
for weakly-ionized atoms, but the collision cross-sections are well-known
for neutral atoms. :numref:`ELSEPAtest` shows a comparison between the model
implemented in Smilei and the theoretical cross-sections calculated by the code
ELSEPA [Salvat2005]_.


.. _ELSEPAtest:

.. figure:: /_static/elsepa.png
  :width: 14cm
  
  Average deflection angle vs. time for collisions between electrons and neutral
  atoms. The lines correspond to Smilei simulations with various electron
  beam energy, and the dots are calculated from the theoretical code ELSEPA.



----

.. _CollNuclearReactions:

.. rst-class:: experimental

Nuclear reactions
^^^^^^^^^^^^^^^^^^^^

Nuclear reactions may occur during collisions when requested. The reaction
scheme is largely inspired from [Higginson2019]_.

.. rubric:: 1. Outline of the nuclear reaction process

We take advantage of the
relativistic kinematics calculations of the binary collision scheme
to introduce the nuclear reactions in the COM frame:

* The cross-section :math:`\sigma` (tabulated for some reactions)
  is interpolated, given the kinetic energies.
* The probability for the reaction to occur is calculated.
* This probability is randomly sampled and, if successful:

  * New macro-particles (the reaction products) are created.
  * Their angle is sampled from a tabulated distribution.
  * Their mpmenta are calculated from the conservation of total energy and momentum.
  * Their momenta are boosted back to the simulation frame.
  
* Otherwise: the collision process proceeds as usual.

.. rubric:: 2. Nuclear reaction probability

The probability for the reaction to occur is calculated as
:math:`P=1-\exp(R\, v\, n\, \sigma\, \Delta t)` where *v* is the relative
velocity, *n* is a corrected density (see [Higginson2020]_), and *R* is
a *rate multiplier* (see [Higginson2019]_).

This factor *R* is of great importance for most applications, because
almost no reactions would occur when :math:`R=1`. This factor artificially
increases the number of reactions to ensure enough statistics. The weights
of the products are adjusted accordingly, and the reactants are not destroyed
in the process: we simply decrease their weight by the same amount.

In Smilei, this factor *R* can be forced by the user to some value, but by
default, it is automatically adjusted so that the final number of created particles
approches the initial number of pairs.

.. rubric:: 3. Creation of the reaction products

Special care must be taken when creating new charged particles while
conserving Poisson's equation. Following Ref. [Higginson2019], we choose to
create two macro-particles of each type. To explain in detail, let us write
the following reaction:

.. math::

  1 + 2 \rightarrow 3 + 4

Two particles of species 3 are created: one at the position of particle 1,
the other at the position of particle 2. Two particles of species 4 are also
created. To conserve the charge at each position, the weights of the new
particles must be:

.. math::

  W_3^{@1} = w \frac{q_1}{q_1+q_2} q_3\\
  W_3^{@2} = w \frac{q_2}{q_1+q_2} q_3\\
  W_4^{@1} = w \frac{q_1}{q_1+q_2} q_4\\
  W_4^{@2} = w \frac{q_2}{q_1+q_2} q_4

where :math:`w` is the products' weight, and the :math:`q_i` are the charges.

.. rubric:: 4. Calculation of the resulting momenta

The conservation of energy reads:

.. math::

  K_1 + K_2 + Q = K_3 + K_4

where the :math:`K_i` are kinetic energies, and :math:`Q` is the reaction's
Q-value. In the COM frame, we have, by definition, equal momenta: :math:`p_3 = p_4`.
Using the relativistic expression :math:`(K_k+m_k)^2=p_k^2+m_k^2`, we can
calculate that

.. math::

  0=p_4^2-p_3^2=K_4 (K_4 + 2m_4) - K_3(K_3+2m_3)

Substituting for :math:`K_4` using the conservation of energy, this translates into

.. math::

  0=A_{00} A_{02} - (A_{20}+A_{02})K_3

where we have defined :math:`A_{ij}=K_1 + K_2 +Q+i\,m_3+j\,m_4`. We thus obtain

.. math::

  K_3 = \frac{A_{00}A_{02}}{A_{20}+A_{02}}\\
  K_3+2m_3 = ... = \frac{A_{20}A_{22}}{A_{20}+A_{02}}

Finally,

.. math:: 

  p_3^2 = K_3(K_3+2m_3) = ... = \frac{A_{00}A_{02}A_{20}A_{22}}{(2A_{11})^2}
  
which expresses the resulting momentum as a function of the initial energies.

----

Collisions debugging
^^^^^^^^^^^^^^^^^^^^

Using the parameter ``debug_every`` in a ``Collisions()`` group (see :ref:`Collisions`)
will create a file with info about these collisions.
These information are stored in the files "Collisions0.h5", "Collisions1.h5", etc.

The *hdf5* files are structured as follows:
  One HDF5 file contains several groups called ``"t********"`` where ``"********"``
  is the timestep. Each of these groups contains several arrays, which represent
  quantities *vs.* space.

The available arrays are:

  * ``s``: defined in [Perez2012]_: :math:`s=N\left<\theta^2\right>`, where :math:`N` is
    the typical number of real collisions during a timestep, and
    :math:`\left<\theta^2\right>` is the average square deviation of individual 
    real collisions. This quantity somewhat represents the typical amount of angular
    deflection accumulated during one timestep.
    **It is recommended that** :math:`s<1` **in order to have realistic collisions.**
  * ``coulomb_log``: average Coulomb logarithm.
  * ``debyelength``: Debye length (not provided if all Coulomb logs are manually defined).

The arrays have the same dimension as the plasma, but each element of these arrays
is an average over all the collisions occurring in a single *patch*.


----

References
^^^^^^^^^^

.. [Desjarlais2001] `M. Desjarlais, Contrib. Plasma Phys. 41, 267 (2001) <http://dx.doi.org/10.1002/1521-3986%28200103%2941%3A2%2F3%3C267%3A%3AAID-CTPP267%3E3.0.CO%3B2-P>`_

.. [Frankel1979] `N. E. Frankel, K. C. Hines, and R. L. Dewar, Phys. Rev. A 20, 2120 (1979) <https://doi.org/10.1103/PhysRevA.20.2120>`_

.. [Higginson2019] `D. P. Higginson, A. Link, A. Schmidt, J. Comput. Phys. 388, 439 (2019) <https://doi.org/10.1016/j.jcp.2019.03.020>`_

.. [Higginson2020] `D. P. Higginson, I. Holod and A. Link, J. Comput. Phys. 413, 109450 (2020) <https://doi.org/10.1016/j.jcp.2020.109450>`_

.. [Lee1984] `Y. T. Lee and R. M. More, Phys. Fluids 27, 1273 (1984) <http://dx.doi.org/10.1063/1.864744>`_

.. [Nanbu1997] `K. Nanbu, Phys. Rev. E 55, 4642 (1997) <http://dx.doi.org/10.1103/PhysRevE.55.4642>`_

.. [Nanbu1998] `K. Nanbu and S. Yonemura, J. Comput. Phys. 145, 639 (1998) <http://dx.doi.org/10.1006/jcph.1998.6049>`_

.. [Perez2012] `F. Pérez et al., Phys. Plasmas 19, 083104 (2012) <http://dx.doi.org/10.1063/1.4742167>`_

.. [Rohrlich1954] `F. Rohrlich and B. C. Carlson, Phys. Rev. 93, 38 (1954) <http://journals.aps.org/pr/abstract/10.1103/PhysRev.93.38>`_

.. [Salvat2005] `F. Salvat, A. Jablonski and C. J. Powell, Comput. Phys. Comm. 165, 157 (2005) <https://doi.org/10.1016/j.cpc.2004.09.006>`_

.. [Sentoku2008] `Y. Sentoku and A. J. Kemp, J. Comput. Phys. 227, 6846 (2008) <http://dx.doi.org/10.1016/j.jcp.2008.03.043>`_


