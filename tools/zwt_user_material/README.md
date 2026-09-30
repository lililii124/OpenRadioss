# ZWT solid user material

Classical Zhu-Wang-Tang nonlinear viscoelasticity through `/MAT/USER01`:
a nonlinear equilibrium spring and two Maxwell branches. Build a double-precision
GNU user library with the official [SDK](https://github.com/OpenRadioss/Tools/tree/main/userlib_sdk)
and load it with `-dylib` in both Starter and Engine. This tool uses the USER01
slot in its own library; it does not assign a built-in LAW number.

## Model and domain

The scalar relation is `sigma = E0*e + alpha*e^2 + beta*e^3 + q1 + q2`,
with `dqi/dt + qi/thetai = Ei*de/dt`; see references [1-2] below,
particularly equations (6)-(7) of [2].
Stress and strain are tensile-positive. Changing from compression-positive
parameters also changes the sign of `alpha`.

The scalar model does not prescribe a unique three-dimensional extension.
Here `F=R U`, `e=U-I`, and `Cnu:e = e/(1+nu) + nu*tr(e)*I/((1+nu)*(1-2nu))`.
In the reference basis, equilibrium stress is `E0*Cnu:e + alpha*e^2 + beta*e^3`
and each branch satisfies `dQi/dt + Qi/thetai = Ei*Cnu:de/dt`. Powers are matrix
products. The returned Cauchy stress is `R*(equilibrium+Q1+Q2)*transpose(R)`.
With `nu=0` and axial strain, this recovers the scalar law exactly.

This is an isotropic **small-stretch** engineering model, allowing superposed
rigid rotations. It requires `||e||F <= emax <= 0.05`, positive `E0` and relaxation
times, nonnegative branch moduli and `-1 < nu < 0.5`. It is not a general
finite-strain thermodynamic model. Damage, temperature coupling, plasticity,
shells and an implicit tangent are not implemented.

Validation requires the sufficient stiffness bound
`E0*min(1/(1+nu),1/(1-2nu)) + min[-emax,emax](2*alpha*x+3*beta*x^2) > 0`,
with a floating-point margin. Some parameter sets outside this certificate
could still be stable; they are deliberately rejected. The wave-speed estimate
uses an upper tangent bound over the entire admitted strain interval.

For strain varying linearly within a step, branches use the exact exponential
update. A series avoids cancellation for small `dt/theta`; the large-ratio
limit avoids overflow (discarded exponential below `exp(-50)`). Zero-time calls
preserve history. Invalid parameters/states preserve input stress and history;
the Engine adapter reports the element/status and stops the unsupported run.

## Input and history

Use isotropic solids with `/PROP/SOLID`: **Ismstr=10, Iframe=1**. The adapter
reads the full global deformation gradient, not incremental engineering shear.
The missing-gradient case is rejected. Other element formulations have not
been validated; the tested configuration is an eight-node one-point brick.

The native density line precedes all **11** constants, in this order:

| Constant | Meaning |
|---|---|
| Kinst, Ginst | `(E0+E1+E2)/(3*(1-2*nu))`, `(E0+E1+E2)/(2*(1+nu))` |
| E0, alpha, beta | Equilibrium coefficients, in stress units |
| E1, theta1 | First branch modulus and relaxation time |
| E2, theta2 | Second branch modulus and relaxation time |
| nu, emax | Common Poisson ratio and maximum strain norm |

For example, synthetic verification parameters in kg, m, s, Pa:

```text
/MAT/USER01/1
ZWT verification material
                1200
7.333333333333333e8 1.1e9 1e9 -5e9 2e10
4e8 1e-2 8e8 1e-4 0
0.05
```

All 18 history variables must be retained across restart: `USR1..6` are `e`,
`USR7..12` are `Q1`, and `USR13..18` are `Q2`, in reference-basis order
`xx, yy, zz, xy, yz, zx`. Strain shear components are tensorial. Standard plastic
strain outputs are zero. A disabled branch (`Ei=0`) must have zero stored stress.
Parameters are verification data, not an experimental material calibration.

## Build and test

From the repository root, with CMake, Ninja and GNU compiler tools on PATH:

```text
git clone https://github.com/OpenRadioss/Tools.git cbuild_zwt/Tools
git -C cbuild_zwt/Tools checkout 4e52942e191d3b1ede4b320fb0f1780f4e41b59a
cmake -S cbuild_zwt/Tools/userlib_sdk/source -B cbuild_zwt/sdk -G Ninja -Darch=win64 -Dcompiler=gfortran -Dprecision=dp
cmake --build cbuild_zwt/sdk
python tools/zwt_user_material/build_userlib.py --sdk-build cbuild_zwt/sdk --build-dir cbuild_zwt/library --verify
```

Use matching GNU compiler versions for SDK modules and the user library.
The helper copies the SDK archive and replaces only its USER01 placeholders.
`--verify` checks analytical ramp/relaxation, extreme time ratios, transactional
rejection and actual Starter/Engine adapters, with floating-point traps enabled.
For Linux use `-Darch=linux64`; that build path is provided but untested.
Windows x64 double precision was checked with GNU Fortran 8.3 and official
OpenRadioss `latest-20260728`. Single precision, MPI and restart-file round trips
have not been validated.

## Model references

1. Wang, L.-L., Huang, D., and Gan, S. (1996). *Nonlinear Viscoelastic
   Constitutive Relations and Nonlinear Viscoelastic Wave Propagation for
   Polymers at High Strain Rates*. In *Constitutive Relation in High/Very High
   Strain Rates*, pp. 137-146. Springer.
   [doi:10.1007/978-4-431-65947-1_16](https://doi.org/10.1007/978-4-431-65947-1_16).
   Model background by Wang and coauthors.
2. Chen, C., Guo, Z., and Tang, E. (2023). *Determination of Elastic Modulus,
   Stress Relaxation Time and Thermal Softening Index in ZWT Constitutive
   Model for Reinforced Al/PTFE*. **Polymers 15**(3), 702.
   [doi:10.3390/polym15030702](https://doi.org/10.3390/polym15030702);
   [open text, section 4.1, equations (6)-(7)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9919274/).
   Only the classical isothermal scalar relation is used; its fitted material
   parameters and thermal extensions are not used here.

The tensor extension, stiffness certificate and numerical integration described
above are implementation choices, not claims about the original scalar model.
