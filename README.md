# ⚗️ CHEMMISOL

CHEMMISOL is a performant chemical equilibrium solver built in the context of the CAMMISOL project.

Widely inspired by the [VMinteq](https://vminteq.com/) software its main purpose is to dynamically solve the equilibrium state of chemical systems.

## Table of contents

- [Feature comparison with VMinteq](#feature-comparison-with-vminteq)
- [Compilation and installation](#compilation-and-installation)
  - [Compile](#compile)
  - [Test](#test)
  - [Install](#install)
- [Specification of chemical systems](#specification-of-chemical-systems)
  - [Reactions](#reactions)
  - [Components and species](#components-and-species)
- [Usage](#usage)
  - [Basic usage](#basic-usage)
  - [Setting and fixing pH and quantities of other chemical components](#setting-and-fixing-ph-and-quantities-of-other-chemical-components)

## Feature comparison with VMinteq

VMinteq is considered as a state of the art reference in the context of chemical speciation computation. However, it is an old fashioned closed source Windows GUI application that cannot be integrated in other models.

The purpose of CHEMMISOL is to provide equivalent features in a flexible, performant, open source, cross platform and embeddable way.

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

# Compilation and installation

## Compile
The library can be compiled as a regular CMAKE project.

```
git clone https://github.com/u2worm/chemmisol-cpp.git
cd chemmisol-cpp
cmake -DCMAKE_BUILD_TYPE=Release -B build -S .
cmake --build build
```

## Test
To run tests:
```
cd build/tests
ctest
```
or:
```
cd build/tests
./chemmisol-tests
```

## Install
The library can be installed on the current system with:
```
cd build
cmake --install .
```

> Depending on your OS and installation, you might need to run `sudo cmake --install .`

# Specification of chemical systems

The definition of chemical systems in CHEMMISOL is based on the concept of components, chemical species and reactions.

Basically, to define a chemical system, it is required to:
- Define a set of reactions.
- Define components and their total quantities.

Components are the building blocks of reactions: the purpose of reactions is to specify how chemical components can combine to form complex chemical species.

## Reactions

Chemical components basically represent canonical species that cannot be divided. A reaction must be specified as a set of reagents, associated to stoichiometric coefficients. By convention:
- Reactants are specified with **positive** coefficients.
- Products are specified with **negative** coefficients.

All reactants and products must correspond to **components** of the chemical system, **except one reagent**, that represents the **produced species** of the reaction.

#### Example 1

Reactions:
1. `H2 <-> H + H`
2. `2 H2O <-> 4 H + O2`

Components:
1. `H`
2. `H2O` (solvent)


The specification of this system is **correct**. The produced species of reaction 1 is then `H2`, and the produced species of reaction 2 is `O2`.

#### Example 2

Reactions:
1. `H2 <-> H + H`
2. `2 H2O <-> 2 H2 + O2`

Components:
1. `H2`
2. `O2`

The specification of this system is **correct**. The produced species of reaction 1 is then `H`, and the produced species of reaction 2 is `H2O`.

#### Example 3

Reactions:
1. `H2 <-> H + H`
2. `2 H2O <-> 4 H + O2`

Components:
1. `H2`
2. `H2O` (solvent)

The specification of this system is **wrong**, since two chemical species that are **not** components are identified in reaction 2: `H` and `O2`. In consequence, the produced species of reaction 2 cannot be properly identified.

## Components and species

A [chemical species](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSpecies.html) is a concrete physical entity that lives in the [chemical system](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html) and can interact with other species. Each chemical species is associated to a quantity (typically in mol), a concentration (typically in mol/l) and an activity (without unit, used in the computation of reaction quotients).

As seen in the reaction example above, chemical species can be defined in two ways:
- Implicitly, as produced species of reactions.
- Explicitly, by declaring components.

Indeed, each [chemical component](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalComponent.html) is automatically associated to a chemical species with the same name. However, the nature of a component is different from the nature of chemical species, since it represents the total quantity of species constituted from this component.

For example, let's consider the following system:

Reactions:
1. `H2O <-> OH- + H+`
2. `PO4 + 3H+ <-> H3PO4`

Components:
1. `H+`
2. `PO4`
3. `H2O` (solvent)

Produced species are then `OH-` and `H3PO4`.

The **total quantity** `N` of the PO4 and H+ components are then defined as follows:

- `N(PO4) = n(PO4) + n(H3PO4)`
- `N(H+) = n(H+) - n(OH-) + 3 * n(H3PO4)`

where n denotes the quantity of each chemical species. Notice that `n(PO4)` and `n(H+)` denote the quantities of PO4 and H+ **species**, that are likely not equal to total quantities of PO4 and H+ **components**. See the [ChemicalComponent](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalComponent.html) documentation for more details about how the total quantities are computed.

In practice, the user input of the model represents the **total quantity** of each component. The actual quantity of each species is then determined by the solved equilibrium state of the chemical system.

# Usage

## Basic usage

The CHEMMISOL API aims to be as simple as possible, in order to facilitate its integration within other projects and models. Here is a basic usage example:

```cpp
#include "chemmisol.h"

using namespace chemmisol;

int main(int, char *[])
{
    ChemicalSystem chemical_system;
    // Defines the reaction H2O <-> OH- + H+ (log K = -13.997)
    chemical_system.addReaction("OH-", -13.997, {
            {"OH-", -1},
            {"H+", -1},
            {"H2O", 1}
            });
    // Defines the reaction Na+ + Cl- <-> NaCl (log K = -0.3)
    chemical_system.addReaction("NaCl", -0.3, {
            {"NaCl", -1},
            {"Na+", AQUEOUS, 1}, // AQUEOUS specified here only for
                                 // demonstration purpose, since it should be
                                 // the default phase
            {"Cl-", AQUEOUS, 1}
            });
    // Defines the reaction H2O + Na+ <-> NaOH + H+ (log K = -13.897)
    chemical_system.addReaction("NaOH", -13.897, {
            {"NaOH", -1},
            {"H+", -1},
            {"Na+", 1},
            {"H2O", 1}
            });

    // Defines the Na+ component and sets its total concentration to 0.1 mol/l
    chemical_system.addComponent("Na+", 0.1*mol/l);
    // Defines the Cl- component and sets its total concentration to 0.1 mol/l
    chemical_system.addComponent("Cl-", 0.1*mol/l);
    // Defines the H2O component as a solvent
    chemical_system.addSolvent("H2O");

    // Automatically adds the H+ component and fixes the pH to 7
    chemical_system.fixPH(7);

    // Solves the equilibrium state
    chemical_system.solveEquilibrium();
    return 0;
}
```

This basic example is available in the `examples` directory and can be run with:
```
cd build
./examples/basic_chemical_system/basic_chemical_system_example
```

Expected output:
```
[chemmisol-core] I Solving chemical equilibrium.
[chemmisol-core] I Init activities:
[chemmisol-core] I  (C) Na+: 0.1
[chemmisol-core] I  (C) Cl-: 0.1
[chemmisol-core] I  (C) H2O: 1
[chemmisol-core] I  (C) H+: 1e-07
[chemmisol-core] I      OH-: 0
[chemmisol-core] I      NaCl: 0
[chemmisol-core] I      NaOH: 0
[chemmisol-core] I Solved activities:
[chemmisol-core] I  (C) Na+: 0.0954352
[chemmisol-core] I  (C) Cl-: 0.0954352
[chemmisol-core] I  (C) H2O: 1
[chemmisol-core] I  (C) H+: 1e-07
[chemmisol-core] I      OH-: 1.00693e-07
[chemmisol-core] I      NaCl: 0.00456476
[chemmisol-core] I      NaOH: 1.20979e-08
```

## Setting and fixing pH and quantities of other chemical components

The pH of the chemical system can be fixed with the [`fixPH()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#ad8f0e5bc375592f1f26a8ccfe15a7249) method.

However, it is also possible to initialize the pH of the system using the [`initPH()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#a8dfa55cc06acd347e941361e40f238b0) method. In this case, the pH is set as the **total quantity** of the `H+` component (or any other user defined component eventually specified as argument by the user). The pH is then dynamically determined depending on the solved chemical equilibrium.

As fixing the pH is equivalent to fixing the concentration of `H+`, the concentration of any component can be fixed in CHEMMISOL, using the [`fixComponent()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#aaa36c3cdf680d4f5f1bae7f1b34eb52e) method.

The total quantity of a component can also be initially specified using the [`addComponent()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#a5a34cbd6f1ed0bbe60ced6a54fb4f919) method.

One of the main advantage of CHEMMISOL is that the total quantities of component, fixed or not, including the pH, can be specified **dynamically**, calling again one of the [`fixPH()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#ad8f0e5bc375592f1f26a8ccfe15a7249), [`fixComponent()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#aaa36c3cdf680d4f5f1bae7f1b34eb52e), [`initPH()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#a8dfa55cc06acd347e941361e40f238b0), [`setTotalQuantity()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#a5a34cbd6f1ed0bbe60ced6a54fb4f919) or [`setTotalConcentration()`](https://u2worm.github.io/chemmisol-cpp/classchemmisol_1_1ChemicalSystem.html#adb72155ff4b60f1d6962f9f9218dc8c9) methods.

```cpp
#include "chemmisol.h"

using namespace chemmisol;

int main(int, char *[])
{
    ChemicalSystem chemical_system;
    // Defines the reaction H2O <-> OH- + H+ (log K = -13.997)
    chemical_system.addReaction("OH-", -13.997, {
            {"OH-", -1},
            {"H+", -1},
            {"H2O", 1}
            });
    // Defines the reaction Na+ + Cl- <-> NaCl (log K = -0.3)
    chemical_system.addReaction("NaCl", -0.3, {
            {"NaCl", -1},
            {"Na+", AQUEOUS, 1}, // AQUEOUS specified here only for
                                 // demonstration purpose, since it should be
                                 // the default phase
            {"Cl-", AQUEOUS, 1}
            });
    // Defines the reaction H2O + Na+ <-> NaOH + H+ (log K = -13.897)
    chemical_system.addReaction("NaOH", -13.897, {
            {"NaOH", -1},
            {"H+", -1},
            {"Na+", 1},
            {"H2O", 1}
            });

    // Defines the Na+ component and sets its total concentration to 0.1 mol/l
    chemical_system.addComponent("Na+", 0.1*mol/l);
    // Defines the Cl- component and sets its total concentration to 0.1 mol/l
    chemical_system.addComponent("Cl-", 0.1*mol/l);
    // Defines the H2O component as a solvent
    chemical_system.addSolvent("H2O");

    // Automatically adds the H+ component and fixes the pH to 7
    chemical_system.fixPH(7);

    // Sets up a basic CSV output
    std::ofstream csv_file("output.csv");
    csv_file << "i, pH, Na+, Cl-, NaCl, NaOH"
        << std::endl;

    for(std::size_t i = 0; i < 10; i++) {
        // Makes the ph vary from 7 to 8,
        chemical_system.fixPH(7 + i*0.1);
        // Reduces the total quantity of Na+ from 0.1 mol/l to 0.09 mol/l
        chemical_system.setTotalConcentration(
                chemical_system.getComponent("Na+"),
                0.1*mol/l - i*0.01*mol/l);

        // Solves the new equilibrium state, from the previous equilibrium
        chemical_system.solveEquilibrium();

        // Concentration output (notice getSpecies() is used, not getComponent())
        csv_file << i << "," <<
            chemical_system.getPH() << "," <<
            chemical_system.getSpecies("Na+").concentration() << "," <<
            chemical_system.getSpecies("Cl-").concentration() << "," <<
            chemical_system.getSpecies("NaCl").concentration() << "," <<
            chemical_system.getSpecies("NaOH").concentration() << std::endl;
    }
}
```

This example is available in the `examples` directory and can be run with:
```
cd build
./examples/dynamic_equilibrium/dynamic_equilibrium_example
```

