Binary collisions
-----------------


Relativistic binary collisions between particles have been implemented in
:program:`Smilei` with the same scheme as the one developed for the
code :program:`Calder`. The following references describe the physics
and numerics of this implementation.

| [Perez2012]_ gives an overview of the technique.
| [Nanbu1997]_ and [Nanbu1998]_ give the original technique from which [Perez2012]_ was developed.
| [Sentoku2008]_, [Lee1984]_ and [Frankel1979]_ provide additional information.

Please refer to :ref:`that doc <Collisions>` for an explanation of how to add collisions in the namelist file.


----

The binary collision scheme
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Collisions are calculated at each timestep and for each collision block
given in the input file.

If *intra-collisions*:
  
  | Create one array of indices pointing to all particles of the species group.
  | Shuffle the array.
  | Split the array in two halves.

If *inter-collisions*:
  
  | Create two arrays of indices pointing to all particles of each species group.
  | Shuffle the largest array. The other array is not shuffled.
  | => The two resulting arrays represent pairs of particles (see algorithm in [Nanbu1998]_).

Calculate a few intermediate quantities:
  
  | Particle density :math:`n_1` of group 1.
  | Particle density :math:`n_2`  of group 2.
  | *Crossed* particle density :math:`n_{12}` (see [Perez2012]_).
  | Other constants.

For each pair of particles:

  | Calculate the momenta in the center-of-mass (COM) frame.
  | Calculate the coulomb log if requested (see [Perez2012]_).
  | Calculate the parameter :math:`s` and its correction at low temperature (see [Perez2012]_).
  | Pick the deflection angle (see [Nanbu1997]_).
  | Deflect particles in the COM frame and go back to the laboratory frame.


.. rubric:: Modification in Smilei

A typo has been introduced in [Perez2012]_, Eq. (22), corresponding to the calculation of
the Coulomb Logarithm: the last parenthesis is written as a squared expression,
but should not. This was corrected in :program:`Smilei`.


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

.. figure:: _static/beam_relaxation123.png
  :width: 10cm
  
  Relaxation of an electron beam. Initial velocity = 0.05, ion charge = 1.
  
.. _beam2:

.. figure:: _static/beam_relaxation456.png
  :width: 10cm
  
  Relaxation of an electron beam. Initial velocity = 0.01, ion charge = 1.

.. _beam3:

.. figure:: _static/beam_relaxation789.png
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

.. figure:: _static/thermalisation_ei123.png
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

.. figure:: _static/temperature_isotropization1.png
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

.. figure:: _static/Maxwellianization1.png
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

.. figure:: _static/Stopping_power123.png
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

.. figure:: _static/conductivity.png
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
  \quad\mathrm{if}\quad 0<k<k_\mathrm{max}
  \\
  \sum\limits_{p=0}^{k-1} \left[ 1+R^{i+k}_{i+p}\left(\frac{W_{i+k}}{W_{i+p}}\bar{P}^{i+p} - \bar{P}^{i+k}\right) \right]
  \prod\limits_{j=0,j\ne p}^{k-1} R^{i+p}_{i+j}
  &
  \quad\mathrm{if}\quad k=k_\mathrm{max}
  \end{array}
  \right.

where :math:`k_\mathrm{max} = Z-Z^\star`.

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

.. figure:: _static/ionization_rate.png
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

.. figure:: _static/ionization_stopping_power.png
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

.. figure:: _static/ionization_multiple.png
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

.. figure:: _static/ionization_recombination.png
  :width: 12cm
  
  Final charge state of various plasmas at various temperatures.

The model does not account for detailed ionization potentials. It provides a rough
approximation, and is particularly questionable for low temperatures or high Z.
We observe that Smilei's approach for impact ionization provides decent estimates
of the ionization state. Detailed comparison to atomic codes has not been done yet.


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

.. [Frankel1979] `N. E. Frankel, K. C. Hines, and R. L. Dewar, Phys. Rev. A 20, 2120 (1979) <http://dx.doi.org/10.1143/JPSJ.67.4084>`_

.. [Lee1984] `Y. T. Lee and R. M. More, Phys. Fluids 27, 1273 (1984) <http://dx.doi.org/10.1063/1.864744>`_

.. [Nanbu1997] `K. Nanbu, Phys. Rev. E 55, 4642 (1997) <http://dx.doi.org/10.1103/PhysRevE.55.4642>`_

.. [Nanbu1998] `K. Nanbu and S. Yonemura, J. Comput. Phys. 145, 639 (1998) <http://dx.doi.org/10.1006/jcph.1998.6049>`_

.. [Perez2012] `F. Pérez et al., Phys. Plasmas 19, 083104 (2012) <http://dx.doi.org/10.1063/1.4742167>`_

.. [Rohrlich1954] `F. Rohrlich and B. C. Carlson, Phys. Rev. 93, 38 (1954) <http://journals.aps.org/pr/abstract/10.1103/PhysRev.93.38>`_

.. [Sentoku2008] `Y. Sentoku and A. J. Kemp, J. Comput. Phys. 227, 6846 (2008) <http://dx.doi.org/10.1016/j.jcp.2008.03.043>`_


