# ⚗️ CHEMMISOL

CHEMMISOL is a performant chemical equilibrium solver built in the context of the CAMMISOL project.

Widely inspired by the [VMinteq](https://vminteq.com/) software its main purpose is to solve the equilibrium state of chemical systems.

## Feature comparison with VMinteq

VMinteq is considered as a state of the art reference in the context of chemical speciation computation. However, it is a closed source Windows GUI application that cannot be integrated in other models.

The purpose of CHEMMISOL is to provide equivalent features in a performant, open source, cross platform and embeddable way.

The following table represents a non-exhaustive list of features supported by VMinteq (except for the "Solver algorithm"), and how they are currently supported by CHEMMISOL.

Even if CHEMMISOL already allows to solve equilibriums in complex chemical systems, a lot of features are still missing, and will probably never be included, depending on our own needs.

<table>
  <tr><th colspan="2"> Feature</th><th>Supported</th><th>Notes</th></tr>
  <tr><td rowspan="2"> Solver algorithm </td><td> Absolute <a href=https://en.wikipedia.org/wiki/Newton%27s_method>Newton-Raphson</a> </td><td> ✔️  </td><td> Even if the equilibrium solving algorithm of VMinteq seems to be based on the Newton method, the exact algorithm used is unclear. In consequence, the "solver algorithm" features cannot be compared directly with VMinteq. </td></tr>
  <tr>                                       <td> <a href=https://en.wikipedia.org/wiki/Numerical_algebraic_geometry>Homotopy continuation</a> + absolute Newton-Raphson </td><td> ❌  </td><td> See [1]. </td></tr>
  <tr><td rowspan="5"> Aqueous system equilibrium solving</td><td>Fixed volume   </td><td> ✔️  </td><td>Equilibrium computed from the <a href=https://en.wikipedia.org/wiki/Conservation_of_mass>law of conservation of mass</a> and equilibrium equations.</td></tr>
  <tr>                                                        <td>Variable volume</td><td> ❌ </td><td> Useful for titration reactions for example. </td></tr>
  <tr>                                                        <td>Infinite solid</td><td> ❌ </td><td> Dissolution of a solid that is never completly dissolved. </td></tr>
  <tr>                                                        <td>Finite solid</td><td> ❌ </td><td> Dissolution of a solid that can be completly dissolved. </td></tr>
  <tr>                                                        <td>Possible solid</td><td> ❌ </td><td> Precipitation of a solid. </td></tr>
  <tr><td rowspan="3"> Activity correction </td><td><a href=https://en.wikipedia.org/wiki/Davies_equation>Davies equations<a></td><td>❌</td><td></td></tr>
  <tr>                                     <td><a href=https://en.wikipedia.org/wiki/Debye%E2%80%93H%C3%BCckel_equation>extended Debye-Hückel equation<a></td><td>❌</td><td></td></tr>
  <tr>                                     <td><a href=https://en.wikipedia.org/wiki/Specific_ion_interaction_theory>SIT<a></td><td>❌</td><td></td></tr>
  <tr><td>Gaseous system</td><td></td><td>❌</td><td></td></tr>
  <tr><td rowspan="4">Mineral adsorption</td><td>Non-electrostatic model (NEM)</td><td>✔️</td><td> Can be used considering the surface complexation model as a pure aqueous chemical system. </td></tr>
  <tr>                                       <td>Diffuse layer model (DLM)</td><td>❌</td><td rowspan="3">See the <a href=https://bit.ly/41NBKEY>VMinteq user guide</a> for detailed explanations.</td></tr>
  <tr>                                       <td>Constant capacitance model (CCM)</td><td>❌</td></tr>
  <tr>                                       <td>Three plane model (TPM)</td><td>❌</td></tr>
</table>

[1] K. Meintjes and A. P. Morgan, “A methodology for solving chemical equilibrium systems,” Applied Mathematics and Computation, vol. 22, no. 4, pp. 333–361, Jun. 1987, doi:[10.1016/0096-3003(87)90076-2](https://doi.org/10.1016/0096-3003(87)90076-2).

