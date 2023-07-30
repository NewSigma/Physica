/*
 * Copyright 2022 WeiBo He.
 *
 * This file is part of Physica.
 *
 * Physica is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Physica is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Physica.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <fstream>
#include "Physica/Core/Interface/VaspWarpper.h"

using namespace Physica::Core;
using ScalarType = Scalar<Float>;
extern const char* outcar;

namespace Physica::Core {
    struct Test {
        static VaspWarpper init() {
            {
                std::ofstream fout("./OUTCAR");
                fout << outcar;
            }

            typename Poscar::LatticeMatrix lattice{3.063970000, 0.000000000, 0.000000000,
                                                  1.355640000, 5.017130000, 0.000000000,
                                                  0.305250000,-1.140470000, 5.195960000};
            typename Poscar::PositionMatrix pos{0.019030000, 0.733400000, 0.307890000,
                                                0.885690000, 0.173690000, 0.053420000,
                                                0.637770000, 0.665250000, 0.810840000,
                                                0.325570000, 0.191560000, 0.558340000,
                                                0.042760000, 0.898360000, 0.689770000,
                                                0.297330000, 0.406950000, 0.941480000,
                                                0.734410000, 0.424230000, 0.428290000,
                                                0.456810000, 0.939370000, 0.176470000};
            VaspWarpper vasp{};
            vasp.vaspWorkingDir = ".";
            vasp.poscar = Poscar({std::move(lattice), std::move(pos), CrystalCell::Type::Direct}, {20, 8}, {2, 6});
            vasp.future = Test::makeDummyFuture();
            return vasp;
        }

        static ProcessFuture makeDummyFuture() {
            ProcessFuture future{};
            future.finished = future.isValid = true;
            return future;
        }
    };
}



int main() {
    VaspWarpper vasp = Physica::Core::Test::init();
    const auto outcar = vasp.getOutcar();
    const auto& force = outcar.getForce();
    const Vector<ScalarType> answer{1.168281, -4.245817, -0.080442, -0.735246, 1.188434, 0.171096, 0.238719, 0.022373, -0.258197, 0.595000, 0.266135, 0.121868, -0.209437, 0.968435, -0.291308, -0.990043, 1.383871, 0.212539, -0.192415, 0.419208, 0.024764, 0.125140, -0.002640, 0.099680};
    if (!vectorNear(force, answer, std::numeric_limits<ScalarType>::epsilon()))
        return 1;
    if (!scalarNear(outcar.getInternalEnergy(), Scalar<Double>(-43.76513486), 1E-8))
      return 1;
    return 0;
}

const char* outcar = " vasp.6.1.0 28Jan20 (build Oct 14 2020 17:24:59) complex                         \n\
   \n\
 executed on             LinuxIFC date 2021.03.18  21:31:16 \n\
 running on   16 total cores \n\
 distrk:  each k-point on   16 cores,    1 groups \n\
 distr:  one band on NCORE=   1 cores,   16 groups \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 INCAR: \n\
 POTCAR:    PAW_PBE Ca_sv 06Sep2000                \n\
 POTCAR:    PAW_PBE O 08Apr2002                    \n\
 ----------------------------------------------------------------------------- \n\
|                                                                             | \n\
|           W    W    AA    RRRRR   N    N  II  N    N   GGGG   !!!           | \n\
|           W    W   A  A   R    R  NN   N  II  NN   N  G    G  !!!           | \n\
|           W    W  A    A  R    R  N N  N  II  N N  N  G       !!!           | \n\
|           W WW W  AAAAAA  RRRRR   N  N N  II  N  N N  G  GGG   !            | \n\
|           WW  WW  A    A  R   R   N   NN  II  N   NN  G    G                | \n\
|           W    W  A    A  R    R  N    N  II  N    N   GGGG   !!!           | \n\
|                                                                             | \n\
|     For optimal performance we recommend to set                             | \n\
|       NCORE = 4 - approx SQRT(number of cores).                             | \n\
|     NCORE specifies how many cores store one orbital (NPAR=cpu/NCORE).      | \n\
|     This setting can greatly improve the performance of VASP for DFT.       | \n\
|     The default, NCORE=1 might be grossly inefficient on modern             | \n\
|     multi-core architectures or massively parallel machines. Do your        | \n\
|     own testing!!!!                                                         | \n\
|     Unfortunately you need to use the default for GW and RPA                | \n\
|     calculations (for HF NCORE is supported but not extensively tested      | \n\
|     yet).                                                                   | \n\
|                                                                             | \n\
 ----------------------------------------------------------------------------- \n\
 \n\
 POTCAR:    PAW_PBE Ca_sv 06Sep2000                \n\
   VRHFIN =Ca: 3s3p4s                                                                                \n\
   LEXCH  = PE                                                                                       \n\
   EATOM  =  1006.0909 eV,   73.9456 Ry                                                              \n\
                                                                                                     \n\
   TITEL  = PAW_PBE Ca_sv 06Sep2000                                                                  \n\
   LULTRA =        F    use ultrasoft PP ?                                                           \n\
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no                                                \n\
   RPACOR =    2.000    partial core radius                                                          \n\
   POMASS =   40.078; ZVAL   =   10.000    mass and valenz                                           \n\
   RCORE  =    2.300    outmost cutoff radius                                                        \n\
   RWIGS  =    2.500; RWIGS  =    1.323    wigner-seitz radius (au A)                                \n\
   ENMAX  =  266.622; ENMIN  =  199.967 eV                                                           \n\
   RCLOC  =    1.808    cutoff for local pot                                                         \n\
   LCOR   =        T    correct aug charges                                                          \n\
   LPAW   =        T    paw PP                                                                       \n\
   EAUG   =  420.852                                                                                 \n\
   RMAX   =    2.359    core radius for proj-oper                                                    \n\
   RAUG   =    1.300    factor for augmentation sphere                                               \n\
   RDEP   =    2.392    radius for radial grids                                                      \n\
   RDEPT  =    1.987    core radius for aug-charge                                                   \n\
                                                                                                     \n\
   Atomic configuration                                                                              \n\
    7 entries                                                                                        \n\
     n  l   j            E        occ.                                                               \n\
     1  0  0.50     -3949.1705   2.0000                                                              \n\
     2  0  0.50      -414.3434   2.0000                                                              \n\
     2  1  1.50      -334.7047   6.0000                                                              \n\
     3  0  0.50       -47.0896   2.0000                                                              \n\
     4  0  0.50        -3.7663   2.0000                                                              \n\
     3  1  1.50       -28.0061   6.0000                                                              \n\
     3  2  1.50        -4.0817   0.0000                                                              \n\
   Description                                                                                       \n\
     l       E           TYP  RCUT    TYP  RCUT                                                      \n\
     0    -47.0896403     23  2.300                                                                  \n\
     0     -3.7663320     23  2.300                                                                  \n\
     1    -28.0060855     23  2.300                                                                  \n\
     1      6.8029130     23  2.300                                                                  \n\
     2     -4.0817478     23  2.300                                                                  \n\
     2      0.4464412     23  2.300                                                                  \n\
  local pseudopotential read in \n\
  partial core-charges read in \n\
  partial kinetic energy density read in \n\
  atomic valenz-charges read in \n\
  non local Contribution for L=           0  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           0  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           1  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           1  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           2  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           2  read in \n\
    real space projection operators read in \n\
    PAW grid and wavefunctions read in \n\
  \n\
   number of l-projection  operators is LMAX  =           6 \n\
   number of lm-projection operators is LMMAX =          18 \n\
  \n\
 POTCAR:    PAW_PBE O 08Apr2002                    \n\
   VRHFIN =O: s2p4                                                                                   \n\
   LEXCH  = PE                                                                                       \n\
   EATOM  =   432.3788 eV,   31.7789 Ry                                                              \n\
                                                                                                     \n\
   TITEL  = PAW_PBE O 08Apr2002                                                                      \n\
   LULTRA =        F    use ultrasoft PP ?                                                           \n\
   IUNSCR =        1    unscreen: 0-lin 1-nonlin 2-no                                                \n\
   RPACOR =    1.200    partial core radius                                                          \n\
   POMASS =   16.000; ZVAL   =    6.000    mass and valenz                                           \n\
   RCORE  =    1.520    outmost cutoff radius                                                        \n\
   RWIGS  =    1.550; RWIGS  =    0.820    wigner-seitz radius (au A)                                \n\
   ENMAX  =  400.000; ENMIN  =  300.000 eV                                                           \n\
   ICORE  =        2    local potential                                                              \n\
   LCOR   =        T    correct aug charges                                                          \n\
   LPAW   =        T    paw PP                                                                       \n\
   EAUG   =  605.392                                                                                 \n\
   DEXC   =    0.000                                                                                 \n\
   RMAX   =    1.553    core radius for proj-oper                                                    \n\
   RAUG   =    1.300    factor for augmentation sphere                                               \n\
   RDEP   =    1.550    radius for radial grids                                                      \n\
   RDEPT  =    1.329    core radius for aug-charge                                                   \n\
                                                                                                     \n\
   Atomic configuration                                                                              \n\
    4 entries                                                                                        \n\
     n  l   j            E        occ.                                                               \n\
     1  0  0.50      -514.6923   2.0000                                                              \n\
     2  0  0.50       -23.9615   2.0000                                                              \n\
     2  1  0.50        -9.0305   4.0000                                                              \n\
     3  2  1.50        -9.5241   0.0000                                                              \n\
   Description                                                                                       \n\
     l       E           TYP  RCUT    TYP  RCUT                                                      \n\
     0    -23.9615318     23  1.200                                                                  \n\
     0     -9.5240782     23  1.200                                                                  \n\
     1     -9.0304911     23  1.520                                                                  \n\
     1      8.1634956     23  1.520                                                                  \n\
     2     -9.5240782      7  1.500                                                                  \n\
  local pseudopotential read in \n\
  partial core-charges read in \n\
  partial kinetic energy density read in \n\
  kinetic energy density of atom read in \n\
  atomic valenz-charges read in \n\
  non local Contribution for L=           0  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           0  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           1  read in \n\
    real space projection operators read in \n\
  non local Contribution for L=           1  read in \n\
    real space projection operators read in \n\
    PAW grid and wavefunctions read in \n\
  \n\
   number of l-projection  operators is LMAX  =           4 \n\
   number of lm-projection operators is LMMAX =           8 \n\
  \n\
  PAW_PBE Ca_sv 06Sep2000               : \n\
 energy of atom  1       EATOM=-1006.0909 \n\
 kinetic energy error for atom=    0.0026 (will be added to EATOM!!) \n\
  PAW_PBE O 08Apr2002                   : \n\
 energy of atom  2       EATOM= -432.3788 \n\
 kinetic energy error for atom=    0.0224 (will be added to EATOM!!) \n\
  \n\
  \n\
 POSCAR: 100000007 \n\
  positions in direct lattice \n\
  No initial velocities read in \n\
 exchange correlation table for  LEXCH =        8 \n\
   RHO(1)=    0.500       N(1)  =     2000 \n\
   RHO(2)=  100.500       N(2)  =     4000 \n\
  \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 ion  position               nearest neighbor table \n\
   1  0.019  0.733  0.308-   8 2.02   5 2.06   8 2.09   7 2.19   6 2.28   7 2.55   4 2.79   2 2.83 \n\
                             4 2.90   3 2.95   3 2.98   3 3.07   2 3.10 \n\
   2  0.886  0.174  0.053-   8 2.07   6 2.08   6 2.10   7 2.12   5 2.12   8 2.16   1 2.83   4 2.92 \n\
                             3 3.01   3 3.02   3 3.02   4 3.07   4 3.08   1 3.10 \n\
   3  0.638  0.665  0.811-   6 2.09   5 2.10   5 2.12   8 2.13   7 2.14   6 2.34   4 2.83   1 2.95 \n\
                             4 2.98   1 2.98   2 3.01   2 3.02   2 3.02   1 3.07 \n\
   4  0.326  0.192  0.558-   6 2.12   7 2.13   7 2.13   5 2.14   8 2.15   5 2.55   1 2.79   3 2.83 \n\
                             1 2.90   2 2.92   3 2.98   2 3.07   2 3.08 \n\
   5  0.043  0.898  0.690-   1 2.06   3 2.10   3 2.12   2 2.12   4 2.14   4 2.55 \n\
   6  0.297  0.407  0.941-   2 2.08   3 2.09   2 2.10   4 2.12   1 2.28   3 2.34 \n\
   7  0.734  0.424  0.428-   2 2.12   4 2.13   4 2.13   3 2.14   1 2.19   1 2.55 \n\
   8  0.457  0.939  0.176-   1 2.02   2 2.07   1 2.09   3 2.13   4 2.15   2 2.16 \n\
  \n\
  LATTYP: Found a triclinic cell. \n\
 ALAT       =     5.2999985685 \n\
 B/A-ratio  =     1.9434860646 \n\
 C/A-ratio  =     0.5781077033 \n\
 COS(alpha) =     0.0046056136 \n\
 COS(beta)  =     0.3223265022 \n\
 COS(gamma) =     0.8188368120 \n\
   \n\
  Lattice vectors: \n\
   \n\
 A1 = (   1.7083300000,  -5.0171300000,   0.0000000000) \n\
 A2 = (   0.0474400000,  -8.8937900000,  -5.1959600000) \n\
 A3 = (   3.0639700000,   0.0000000000,   0.0000000000) \n\
 \n\
 \n\
Analysis of symmetry for initial positions (statically): \n\
===================================================================== \n\
 Subroutine PRICEL returns: \n\
 Original cell was already a primitive cell. \n\
  \n\
 \n\
 Routine SETGRP: Setting up the symmetry group for a  \n\
 triclinic supercell. \n\
 \n\
 \n\
 Subroutine GETGRP returns: Found  1 space group operations \n\
 (whereof  1 operations were pure point group operations) \n\
 out of a pool of  2 trial point group operations. \n\
 \n\
 \n\
The static configuration has the point symmetry C_1 . \n\
 \n\
 \n\
Analysis of symmetry for dynamics (positions and initial velocities): \n\
===================================================================== \n\
 Subroutine PRICEL returns: \n\
 Original cell was already a primitive cell. \n\
  \n\
 \n\
 Routine SETGRP: Setting up the symmetry group for a  \n\
 triclinic supercell. \n\
 \n\
 \n\
 Subroutine GETGRP returns: Found  1 space group operations \n\
 (whereof  1 operations were pure point group operations) \n\
 out of a pool of  2 trial point group operations. \n\
 \n\
 \n\
The dynamic configuration has the point symmetry C_1 . \n\
 \n\
 \n\
 Subroutine INISYM returns: Found  1 space group operations \n\
 (whereof  1 operations are pure point group operations), \n\
 and found     1 'primitive' translations \n\
 \n\
 \n\
---------------------------------------------------------------------------------------- \n\
 \n\
                                     Primitive cell                                      \n\
 \n\
  volume of cell :      79.8740 \n\
 \n\
  direct lattice vectors                    reciprocal lattice vectors \n\
     3.063970000  0.000000000  0.000000000     0.326373953 -0.088186988 -0.038529986 \n\
     1.355640000  5.017130000  0.000000000     0.000000000  0.199317139  0.043748454 \n\
     0.305250000 -1.140470000  5.195960000     0.000000000  0.000000000  0.192457217 \n\
 \n\
  length of vectors \n\
     3.063970000  5.197052361  5.328400295     0.340266751  0.204061876  0.192457217 \n\
 \n\
  position of ions in fractional coordinates (direct lattice) \n\
     0.019030000  0.733400000  0.307890000 \n\
     0.885690000  0.173690000  0.053420000 \n\
     0.637770000  0.665250000  0.810840000 \n\
     0.325570000  0.191560000  0.558340000 \n\
     0.042760000  0.898360000  0.689770000 \n\
     0.297330000  0.406950000  0.941480000 \n\
     0.734410000  0.424230000  0.428290000 \n\
     0.456810000  0.939370000  0.176470000 \n\
 \n\
  ion indices of the primitive-cell ions \n\
   primitive index   ion index \n\
                 1           1 \n\
                 2           2 \n\
                 3           3 \n\
                 4           4 \n\
                 5           5 \n\
                 6           6 \n\
                 7           7 \n\
                 8           8 \n\
 \n\
---------------------------------------------------------------------------------------- \n\
 \n\
  \n\
  \n\
 KPOINTS: Monkhorst-Pack                           \n\
 \n\
Automatic generation of k-mesh. \n\
 generate k-points for:    3    2    2 \n\
Space group operators: \n\
 irot       det(A)        alpha          n_x          n_y          n_z        tau_x        tau_y        tau_z \n\
    1     1.000000     0.000000     1.000000     0.000000     0.000000     0.000000     0.000000     0.000000 \n\
  \n\
 Subroutine IBZKPT returns following result: \n\
 =========================================== \n\
  \n\
 Found      8 irreducible k-points: \n\
  \n\
 Following reciprocal coordinates: \n\
            Coordinates               Weight \n\
  0.000000  0.000000  0.000000      1.000000 \n\
  0.333333  0.000000  0.000000      2.000000 \n\
  0.000000  0.500000 -0.000000      1.000000 \n\
  0.333333  0.500000 -0.000000      2.000000 \n\
  0.000000  0.000000  0.500000      1.000000 \n\
  0.333333  0.000000  0.500000      2.000000 \n\
  0.000000  0.500000  0.500000      1.000000 \n\
  0.333333  0.500000  0.500000      2.000000 \n\
  \n\
 Following cartesian coordinates: \n\
            Coordinates               Weight \n\
  0.000000  0.000000  0.000000      1.000000 \n\
  0.108791 -0.029396 -0.012843      2.000000 \n\
  0.000000  0.099659  0.021874      1.000000 \n\
  0.108791  0.070263  0.009031      2.000000 \n\
  0.000000  0.000000  0.096229      1.000000 \n\
  0.108791 -0.029396  0.083385      2.000000 \n\
  0.000000  0.099659  0.118103      1.000000 \n\
  0.108791  0.070263  0.105260      2.000000 \n\
  \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
 Dimension of arrays: \n\
   k-points           NKPTS =      8   k-points in BZ     NKDIM =      8   number of bands    NBANDS=     48 \n\
   number of dos      NEDOS =    301   number of ions     NIONS =      8 \n\
   non local maximal  LDIM  =      6   non local SUM 2l+1 LMDIM =     18 \n\
   total plane-waves  NPLWV =  15120 \n\
   max r-space proj   IRMAX =      1   max aug-charges    IRDMAX=  72445 \n\
   dimension x,y,z NGX =    18 NGY =   28 NGZ =   30 \n\
   dimension x,y,z NGXF=    36 NGYF=   56 NGZF=   60 \n\
   support grid    NGXF=    72 NGYF=  112 NGZF=  120 \n\
   ions per type =               4   4 \n\
   NGX,Y,Z   is equivalent  to a cutoff of   9.77,  8.96,  9.36 a.u. \n\
   NGXF,Y,Z  is equivalent  to a cutoff of  19.53, 17.91, 18.72 a.u. \n\
 \n\
 SYSTEM =  unknown system                           \n\
 POSCAR =  100000007                                \n\
 \n\
 Startparameter for this run: \n\
   NWRITE =      2    write-flag & timer \n\
   PREC   = normal    normal or accurate (medium, high low for compatibility) \n\
   ISTART =      0    job   : 0-new  1-cont  2-samecut \n\
   ICHARG =      2    charge: 1-file 2-atom 10-const \n\
   ISPIN  =      1    spin polarized calculation? \n\
   LNONCOLLINEAR =      F non collinear calculations \n\
   LSORBIT =      F    spin-orbit coupling \n\
   INIWAV =      1    electr: 0-lowe 1-rand  2-diag \n\
   LASPH  =      F    aspherical Exc in radial PAW \n\
   METAGGA=      F    non-selfconsistent MetaGGA calc. \n\
 \n\
 Electronic Relaxation 1 \n\
   ENCUT  =  500.0 eV  36.75 Ry    6.06 a.u.   5.59  9.48  9.71*2*pi/ulx,y,z \n\
   ENINI  =  500.0     initial cutoff \n\
   ENAUG  =  605.4 eV  augmentation charge cutoff \n\
   NELM   =     60;   NELMIN=  2; NELMDL= -5     # of ELM steps  \n\
   EDIFF  = 0.1E-04   stopping-criterion for ELM \n\
   LREAL  =      F    real-space projection \n\
   NLSPLINE    = F    spline interpolate recip. space projectors \n\
   LCOMPAT=      F    compatible to vasp.4.4 \n\
   GGA_COMPAT  = T    GGA compatible to vasp.4.4-vasp.4.6 \n\
   LMAXPAW     = -100 max onsite density \n\
   LMAXMIX     =    2 max onsite mixed and CHGCAR \n\
   VOSKOWN=      0    Vosko Wilk Nusair interpolation \n\
   ROPT   =    0.00000   0.00000 \n\
 Ionic relaxation \n\
   EDIFFG = 0.1E-03   stopping-criterion for IOM \n\
   NSW    =      0    number of steps for IOM \n\
   NBLOCK =      1;   KBLOCK =      1    inner block; outer block  \n\
   IBRION =     -1    ionic relax: 0-MD 1-quasi-New 2-CG \n\
   NFREE  =      0    steps in history (QN), initial steepest desc. (CG) \n\
   ISIF   =      2    stress and relaxation \n\
   IWAVPR =     10    prediction:  0-non 1-charg 2-wave 3-comb \n\
   ISYM   =      2    0-nonsym 1-usesym 2-fastsym \n\
   LCORR  =      T    Harris-Foulkes like correction to forces \n\
 \n\
   POTIM  = 0.5000    time-step for ionic-motion \n\
   TEIN   =    0.0    initial temperature \n\
   TEBEG  =    0.0;   TEEND  =   0.0 temperature during run \n\
   SMASS  =  -3.00    Nose mass-parameter (am) \n\
   estimated Nose-frequenzy (Omega)   =  0.10E-29 period in steps = 0.13E+47 mass=  -0.215E-27a.u. \n\
   SCALEE = 1.0000    scale energy and forces \n\
   NPACO  =    256;   APACO  = 16.0  distance and # of slots for P.C. \n\
   PSTRESS=  650.0 pullay stress \n\
 \n\
  Mass of Ions in am \n\
   POMASS =  40.08 16.00 \n\
  Ionic Valenz \n\
   ZVAL   =  10.00  6.00 \n\
  Atomic Wigner-Seitz radii \n\
   RWIGS  =  -1.00 -1.00 \n\
  virtual crystal weights  \n\
   VCA    =   1.00  1.00 \n\
   NELECT =      64.0000    total number of electrons \n\
   NUPDOWN=      -1.0000    fix difference up-down \n\
 \n\
 DOS related values: \n\
   EMIN   =  10.00;   EMAX   =-10.00  energy-range for DOS \n\
   EFERMI =   0.00 \n\
   ISMEAR =     0;   SIGMA  =   0.05  broadening in eV -4-tet -1-fermi 0-gaus \n\
 \n\
 Electronic relaxation 2 (details) \n\
   IALGO  =     38    algorithm \n\
   LDIAG  =      T    sub-space diagonalisation (order eigenvalues) \n\
   LSUBROT=      F    optimize rotation matrix (better conditioning) \n\
   TURBO    =      0    0=normal 1=particle mesh \n\
   IRESTART =      0    0=no restart 2=restart with 2 vectors \n\
   NREBOOT  =      0    no. of reboots \n\
   NMIN     =      0    reboot dimension \n\
   EREF     =   0.00    reference energy to select bands \n\
   IMIX   =      4    mixing-type and parameters \n\
     AMIX     =   0.40;   BMIX     =  1.00 \n\
     AMIX_MAG =   1.60;   BMIX_MAG =  1.00 \n\
     AMIN     =   0.10 \n\
     WC   =   100.;   INIMIX=   1;  MIXPRE=   1;  MAXMIX= -45 \n\
 \n\
 Intra band minimization: \n\
   WEIMIN = 0.0000     energy-eigenvalue tresh-hold \n\
   EBREAK =  0.52E-07  absolut break condition \n\
   DEPER  =   0.30     relativ break condition   \n\
 \n\
   TIME   =   0.40     timestep for ELM \n\
 \n\
  volume/ion in A,a.u.               =       9.98        67.38 \n\
  Fermi-wavevector in a.u.,A,eV,Ry     =   1.520546  2.873416 31.457494  2.312061 \n\
  Thomas-Fermi vector in A             =   2.629382 \n\
  \n\
 Write flags \n\
   LWAVE        =      T    write WAVECAR \n\
   LDOWNSAMPLE  =      F    k-point downsampling of WAVECAR \n\
   LCHARG       =      T    write CHGCAR \n\
   LVTOT        =      F    write LOCPOT, total local potential \n\
   LVHAR        =      F    write LOCPOT, Hartree potential only \n\
   LELF         =      F    write electronic localiz. function (ELF) \n\
   LORBIT       =      0    0 simple, 1 ext, 2 COOP (PROOUT), +10 PAW based schemes \n\
 \n\
 \n\
 Dipole corrections \n\
   LMONO  =      F    monopole corrections only (constant potential shift) \n\
   LDIPOL =      F    correct potential (dipole corrections) \n\
   IDIPOL =      0    1-x, 2-y, 3-z, 4-all directions  \n\
   EPSILON=  1.0000000 bulk dielectric constant \n\
 \n\
 Exchange correlation treatment: \n\
   GGA     =    --    GGA type \n\
   LEXCH   =     8    internal setting for exchange type \n\
   VOSKOWN=      0    Vosko Wilk Nusair interpolation \n\
   LHFCALC =     F    Hartree Fock is set to \n\
   LHFONE  =     F    Hartree Fock one center treatment \n\
   AEXX    =    0.0000 exact exchange contribution \n\
 \n\
 Linear response parameters \n\
   LEPSILON=     F    determine dielectric tensor \n\
   LRPA    =     F    only Hartree local field effects (RPA) \n\
   LNABLA  =     F    use nabla operator in PAW spheres \n\
   LVEL    =     F    velocity operator in full k-point grid \n\
   LINTERFAST=   F  fast interpolation \n\
   KINTER  =     0    interpolate to denser k-point grid \n\
   CSHIFT  =0.1000    complex shift for real part using Kramers Kronig \n\
   OMEGAMAX=  -1.0    maximum frequency \n\
   DEG_THRESHOLD= 0.2000000E-02 threshold for treating states as degnerate \n\
   RTIME   =   -0.100 relaxation time in fs \n\
  (WPLASMAI=    0.000 imaginary part of plasma frequency in eV, 0.658/RTIME) \n\
   DFIELD  = 0.0000000 0.0000000 0.0000000 field for delta impulse in time \n\
  \n\
 Orbital magnetization related: \n\
   ORBITALMAG=     F  switch on orbital magnetization \n\
   LCHIMAG   =     F  perturbation theory with respect to B field \n\
   DQ        =  0.001000  dq finite difference perturbation B field \n\
   LLRAUG    =     F  two centre corrections for induced B field \n\
 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 Static calculation \n\
 charge density and potential will be updated during run \n\
 non-spin polarized calculation \n\
 Variant of blocked Davidson \n\
 Davidson routine will perform the subspace rotation \n\
 perform sub-space diagonalisation \n\
    after iterative eigenvector-optimisation \n\
 modified Broyden-mixing scheme, WC =      100.0 \n\
 initial mixing is a Kerker type mixing with AMIX =  0.4000 and BMIX =      1.0000 \n\
 Hartree-type preconditioning will be used \n\
 using additional bands           16 \n\
 reciprocal scheme for non local part \n\
 use partial core corrections \n\
 calculate Harris-corrections to forces  \n\
   (improved forces if not selfconsistent) \n\
 use gradient corrections  \n\
 use of overlap-Matrix (Vanderbilt PP) \n\
 Gauss-broadening in eV      SIGMA  =   0.05 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
  energy-cutoff  :      500.00 \n\
  volume of cell :       79.87 \n\
      direct lattice vectors                 reciprocal lattice vectors \n\
     3.063970000  0.000000000  0.000000000     0.326373953 -0.088186988 -0.038529986 \n\
     1.355640000  5.017130000  0.000000000     0.000000000  0.199317139  0.043748454 \n\
     0.305250000 -1.140470000  5.195960000     0.000000000  0.000000000  0.192457217 \n\
 \n\
  length of vectors \n\
     3.063970000  5.197052361  5.328400295     0.340266751  0.204061876  0.192457217 \n\
 \n\
 \n\
  \n\
 k-points in units of 2pi/SCALE and weight: Monkhorst-Pack                           \n\
   0.00000000  0.00000000  0.00000000       0.083 \n\
   0.10879132 -0.02939566 -0.01284333       0.167 \n\
   0.00000000  0.09965857  0.02187423       0.083 \n\
   0.10879132  0.07026291  0.00903090       0.167 \n\
   0.00000000  0.00000000  0.09622861       0.083 \n\
   0.10879132 -0.02939566  0.08338528       0.167 \n\
   0.00000000  0.09965857  0.11810284       0.083 \n\
   0.10879132  0.07026291  0.10525951       0.167 \n\
  \n\
 k-points in reciprocal lattice and weights: Monkhorst-Pack                           \n\
   0.00000000  0.00000000  0.00000000       0.083 \n\
   0.33333333  0.00000000  0.00000000       0.167 \n\
   0.00000000  0.50000000 -0.00000000       0.083 \n\
   0.33333333  0.50000000 -0.00000000       0.167 \n\
   0.00000000  0.00000000  0.50000000       0.083 \n\
   0.33333333  0.00000000  0.50000000       0.167 \n\
   0.00000000  0.50000000  0.50000000       0.083 \n\
   0.33333333  0.50000000  0.50000000       0.167 \n\
  \n\
 position of ions in fractional coordinates (direct lattice)  \n\
   0.01903000  0.73340000  0.30789000 \n\
   0.88569000  0.17369000  0.05342000 \n\
   0.63777000  0.66525000  0.81084000 \n\
   0.32557000  0.19156000  0.55834000 \n\
   0.04276000  0.89836000  0.68977000 \n\
   0.29733000  0.40695000  0.94148000 \n\
   0.73441000  0.42423000  0.42829000 \n\
   0.45681000  0.93937000  0.17647000 \n\
  \n\
 position of ions in cartesian coordinates  (Angst): \n\
   1.14651715  3.32842383  1.59978412 \n\
   2.96549516  0.81050140  0.27756818 \n\
   3.10345657  2.41290704  4.21309221 \n\
   1.42765640  0.32431140  2.90111231 \n\
   1.55942040  3.72052691  3.58401733 \n\
   1.75007467  0.96799136  4.89189242 \n\
   2.95604889  1.63996516  2.22537771 \n\
   2.72696715  4.51168267  0.91693106 \n\
  \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 k-point   1 :   0.0000 0.0000 0.0000  plane waves:    2037 \n\
 k-point   2 :   0.3333 0.0000 0.0000  plane waves:    2024 \n\
 k-point   3 :   0.0000 0.5000-0.0000  plane waves:    2036 \n\
 k-point   4 :   0.3333 0.5000-0.0000  plane waves:    2027 \n\
 k-point   5 :   0.0000 0.0000 0.5000  plane waves:    2030 \n\
 k-point   6 :   0.3333 0.0000 0.5000  plane waves:    2031 \n\
 k-point   7 :   0.0000 0.5000 0.5000  plane waves:    2030 \n\
 k-point   8 :   0.3333 0.5000 0.5000  plane waves:    2015 \n\
 \n\
 maximum and minimum number of plane-waves per node :      2037     2015 \n\
 \n\
 maximum number of plane-waves:      2037 \n\
 maximum index in each direction:  \n\
   IXMAX=    5   IYMAX=    9   IZMAX=    9 \n\
   IXMIN=   -5   IYMIN=   -9   IZMIN=  -10 \n\
 \n\
 \n\
 serial   3D FFT for wavefunctions \n\
 parallel 3D FFT for charge: \n\
    minimum data exchange during FFTs selected (reduces bandwidth) \n\
 \n\
 \n\
 total amount of memory used by VASP MPI-rank0    38421. kBytes \n\
======================================================================= \n\
 \n\
   base      :      30000. kBytes \n\
   nonl-proj :       4980. kBytes \n\
   fftplans  :       1298. kBytes \n\
   grid      :       1190. kBytes \n\
   one-center:        124. kBytes \n\
   wavefun   :        829. kBytes \n\
  \n\
     INWAV:  cpu time    0.0000: real time    0.0002 \n\
 Broyden mixing: mesh for mixing (old mesh) \n\
   NGX = 11   NGY = 19   NGZ = 19 \n\
  (NGX  = 36   NGY  = 56   NGZ  = 60) \n\
  gives a total of   3971 points \n\
 \n\
 initial charge density was supplied: \n\
 charge density of overlapping atoms calculated \n\
 number of electron      64.0000000 magnetization  \n\
 keeping initial charge density in first step \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 Maximum index for augmentation-charges         5360 (set IRDMAX) \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 First call to EWALD:  gamma=   0.412 \n\
 Maximum number of real-space cells 4x 2x 2 \n\
 Maximum number of reciprocal cells 2x 3x 3 \n\
 \n\
    FEWALD:  cpu time    0.0020: real time    0.0024 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   1)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0210: real time    0.0229 \n\
    SETDIJ:  cpu time    0.0760: real time    0.0751 \n\
     EDDAV:  cpu time    0.6659: real time    0.6674 \n\
       DOS:  cpu time    0.0010: real time    0.0011 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.7639: real time    0.7665 \n\
 \n\
 eigenvalue-minimisations  :   768 \n\
 total energy-change (2. order) : 0.6468020E+03  (-0.3255381E+04) \n\
 number of electron      64.0000000 magnetization  \n\
 augmentation part       64.0000000 magnetization  \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -923.94447659 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       290.06244233 \n\
  PAW double counting   =      4370.59195404    -4337.70847039 \n\
  entropy T*S    EENTRO =        -0.00007462 \n\
  eigenvalues    EBANDS =        37.82325063 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       646.80196115 eV \n\
 \n\
  energy without entropy =      646.80203577  energy(sigma->0) =      646.80199846 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   2)  --------------------------------------- \n\
 \n\
 \n\
     EDDAV:  cpu time    0.7779: real time    0.7796 \n\
       DOS:  cpu time    0.0010: real time    0.0004 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.7789: real time    0.7800 \n\
 \n\
 eigenvalue-minimisations  :   912 \n\
 total energy-change (2. order) :-0.6547233E+03  (-0.6172637E+03) \n\
 number of electron      64.0000000 magnetization  \n\
 augmentation part       64.0000000 magnetization  \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -923.94447659 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       290.06244233 \n\
  PAW double counting   =      4370.59195404    -4337.70847039 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -616.90010719 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =        -7.92132205 eV \n\
 \n\
  energy without entropy =       -7.92132205  energy(sigma->0) =       -7.92132205 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   3)  --------------------------------------- \n\
 \n\
 \n\
     EDDAV:  cpu time    1.1318: real time    1.1330 \n\
       DOS:  cpu time    0.0000: real time    0.0004 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    1.1318: real time    1.1334 \n\
 \n\
 eigenvalue-minimisations  :  1280 \n\
 total energy-change (2. order) :-0.4234436E+02  (-0.4204828E+02) \n\
 number of electron      64.0000000 magnetization  \n\
 augmentation part       64.0000000 magnetization  \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -923.94447659 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       290.06244233 \n\
  PAW double counting   =      4370.59195404    -4337.70847039 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -659.24446471 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -50.26567957 eV \n\
 \n\
  energy without entropy =      -50.26567957  energy(sigma->0) =      -50.26567957 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   4)  --------------------------------------- \n\
 \n\
 \n\
     EDDAV:  cpu time    0.7949: real time    0.7954 \n\
       DOS:  cpu time    0.0010: real time    0.0006 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.7959: real time    0.7960 \n\
 \n\
 eigenvalue-minimisations  :   928 \n\
 total energy-change (2. order) :-0.1616046E+01  (-0.1615530E+01) \n\
 number of electron      64.0000000 magnetization  \n\
 augmentation part       64.0000000 magnetization  \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -923.94447659 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       290.06244233 \n\
  PAW double counting   =      4370.59195404    -4337.70847039 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -660.86051102 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -51.88172589 eV \n\
 \n\
  energy without entropy =      -51.88172589  energy(sigma->0) =      -51.88172589 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   5)  --------------------------------------- \n\
 \n\
 \n\
     EDDAV:  cpu time    0.9919: real time    0.9945 \n\
       DOS:  cpu time    0.0010: real time    0.0004 \n\
    CHARGE:  cpu time    0.0700: real time    0.0713 \n\
    MIXING:  cpu time    0.0020: real time    0.0016 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    1.0648: real time    1.0678 \n\
 \n\
 eigenvalue-minimisations  :  1152 \n\
 total energy-change (2. order) :-0.1127044E-01  (-0.1127019E-01) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       13.8419675 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.18455E+01    rms(broyden)= 0.18451E+01 \n\
  rms(prec ) = 0.35170E+01 \n\
  weight for this iteration     100.00 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -923.94447659 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       290.06244233 \n\
  PAW double counting   =      4370.59195404    -4337.70847039 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -660.87178146 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -51.89299632 eV \n\
 \n\
  energy without entropy =      -51.89299632  energy(sigma->0) =      -51.89299632 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   6)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0150: real time    0.0151 \n\
    SETDIJ:  cpu time    0.0460: real time    0.0460 \n\
     EDDAV:  cpu time    1.0598: real time    1.0614 \n\
       DOS:  cpu time    0.0000: real time    0.0004 \n\
    CHARGE:  cpu time    0.0700: real time    0.0704 \n\
    MIXING:  cpu time    0.0020: real time    0.0017 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    1.1928: real time    1.1949 \n\
 \n\
 eigenvalue-minimisations  :  1200 \n\
 total energy-change (2. order) : 0.8223664E+01  (-0.3158256E+01) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       13.3746161 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.58836E+00    rms(broyden)= 0.58746E+00 \n\
  rms(prec ) = 0.79543E+00 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   0.9020 \n\
  0.9020 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -972.97432548 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       293.93888016 \n\
  PAW double counting   =      4695.32559559    -4678.96068405 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -590.97613460 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.66933263 eV \n\
 \n\
  energy without entropy =      -43.66933263  energy(sigma->0) =      -43.66933263 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   7)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0160: real time    0.0158 \n\
    SETDIJ:  cpu time    0.0460: real time    0.0458 \n\
     EDDAV:  cpu time    0.7699: real time    0.7701 \n\
       DOS:  cpu time    0.0000: real time    0.0004 \n\
    CHARGE:  cpu time    0.0530: real time    0.0537 \n\
    MIXING:  cpu time    0.0020: real time    0.0017 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.8869: real time    0.8874 \n\
 \n\
 eigenvalue-minimisations  :   944 \n\
 total energy-change (2. order) :-0.1260917E+00  (-0.4221306E+00) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       13.0344358 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.35967E+00    rms(broyden)= 0.35965E+00 \n\
  rms(prec ) = 0.48779E+00 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   1.3941 \n\
  1.0748  1.7135 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -971.20450891 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       293.81037887 \n\
  PAW double counting   =      4864.56226702    -4855.28596005 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -585.65493705 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.79542438 eV \n\
 \n\
  energy without entropy =      -43.79542438  energy(sigma->0) =      -43.79542438 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   8)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0140: real time    0.0137 \n\
    SETDIJ:  cpu time    0.0320: real time    0.0321 \n\
     EDDAV:  cpu time    0.8039: real time    0.8044 \n\
       DOS:  cpu time    0.0000: real time    0.0003 \n\
    CHARGE:  cpu time    0.0400: real time    0.0402 \n\
    MIXING:  cpu time    0.0020: real time    0.0013 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.8919: real time    0.8920 \n\
 \n\
 eigenvalue-minimisations  :  1152 \n\
 total energy-change (2. order) : 0.1675737E-01  (-0.4334096E-01) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       13.0243194 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.11334E+00    rms(broyden)= 0.11330E+00 \n\
  rms(prec ) = 0.22505E+00 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   1.4727 \n\
  2.4770  0.9706  0.9706 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -969.62731590 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       293.69881775 \n\
  PAW double counting   =      5139.02343216    -5134.52555016 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -582.32538662 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.77866701 eV \n\
 \n\
  energy without entropy =      -43.77866701  energy(sigma->0) =      -43.77866701 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(   9)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0120: real time    0.0128 \n\
    SETDIJ:  cpu time    0.0240: real time    0.0243 \n\
     EDDAV:  cpu time    0.5219: real time    0.5216 \n\
       DOS:  cpu time    0.0000: real time    0.0004 \n\
    CHARGE:  cpu time    0.0460: real time    0.0467 \n\
    MIXING:  cpu time    0.0020: real time    0.0015 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.6059: real time    0.6074 \n\
 \n\
 eigenvalue-minimisations  :   864 \n\
 total energy-change (2. order) : 0.1572483E-01  (-0.1250375E-01) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       12.9930264 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.21115E-01    rms(broyden)= 0.21043E-01 \n\
  rms(prec ) = 0.31584E-01 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   1.4871 \n\
  2.3439  1.7085  0.9479  0.9479 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -973.33666088 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       293.99705646 \n\
  PAW double counting   =      5234.31918947    -5231.95409428 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -576.76576870 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.76294217 eV \n\
 \n\
  energy without entropy =      -43.76294217  energy(sigma->0) =      -43.76294217 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(  10)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0130: real time    0.0133 \n\
    SETDIJ:  cpu time    0.0280: real time    0.0276 \n\
     EDDAV:  cpu time    0.7799: real time    0.7811 \n\
       DOS:  cpu time    0.0010: real time    0.0004 \n\
    CHARGE:  cpu time    0.0440: real time    0.0446 \n\
    MIXING:  cpu time    0.0020: real time    0.0015 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.8679: real time    0.8684 \n\
 \n\
 eigenvalue-minimisations  :  1104 \n\
 total energy-change (2. order) :-0.1629093E-02  (-0.1388525E-02) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       12.9719221 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.15127E-01    rms(broyden)= 0.15127E-01 \n\
  rms(prec ) = 0.21468E-01 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   1.5263 \n\
  2.4777  2.4777  0.9065  0.9065  0.8630 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -973.28343560 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       293.98754065 \n\
  PAW double counting   =      5225.73341071    -5222.98102914 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -577.19839364 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.76457127 eV \n\
 \n\
  energy without entropy =      -43.76457127  energy(sigma->0) =      -43.76457127 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(  11)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0140: real time    0.0142 \n\
    SETDIJ:  cpu time    0.0270: real time    0.0275 \n\
     EDDAV:  cpu time    0.6699: real time    0.6707 \n\
       DOS:  cpu time    0.0000: real time    0.0004 \n\
    CHARGE:  cpu time    0.0450: real time    0.0446 \n\
    MIXING:  cpu time    0.0020: real time    0.0015 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.7579: real time    0.7589 \n\
 \n\
 eigenvalue-minimisations  :  1008 \n\
 total energy-change (2. order) :-0.4843882E-03  (-0.5550439E-04) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       12.9710026 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.68573E-02    rms(broyden)= 0.68567E-02 \n\
  rms(prec ) = 0.10444E-01 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   1.5344 \n\
  2.7259  2.3135  1.3533  0.9336  0.9336  0.9466 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -973.63548355 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       294.00662105 \n\
  PAW double counting   =      5225.95737115    -5222.81517356 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -577.25572649 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.76505565 eV \n\
 \n\
  energy without entropy =      -43.76505565  energy(sigma->0) =      -43.76505565 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(  12)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0140: real time    0.0142 \n\
    SETDIJ:  cpu time    0.0280: real time    0.0275 \n\
     EDDAV:  cpu time    0.6689: real time    0.6692 \n\
       DOS:  cpu time    0.0000: real time    0.0004 \n\
    CHARGE:  cpu time    0.0450: real time    0.0448 \n\
    MIXING:  cpu time    0.0010: real time    0.0015 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.7569: real time    0.7576 \n\
 \n\
 eigenvalue-minimisations  :   992 \n\
 total energy-change (2. order) :-0.7583995E-04  (-0.1299776E-04) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       12.9708022 magnetization  \n\
 \n\
 Broyden mixing: \n\
  rms(total) = 0.18885E-02    rms(broyden)= 0.18856E-02 \n\
  rms(prec ) = 0.27840E-02 \n\
  weight for this iteration     100.00 \n\
 \n\
 eigenvalues of (default mixing * dielectric matrix) \n\
  average eigenvalue GAMMA=   1.5239 \n\
  2.8564  2.3091  1.7761  0.9220  0.9220  0.9409  0.9409 \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -973.87272525 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       294.02051111 \n\
  PAW double counting   =      5229.75122134    -5226.48931138 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -577.15216307 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.76513149 eV \n\
 \n\
  energy without entropy =      -43.76513149  energy(sigma->0) =      -43.76513149 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
--------------------------------------- Iteration      1(  13)  --------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0140: real time    0.0142 \n\
    SETDIJ:  cpu time    0.0270: real time    0.0277 \n\
     EDDAV:  cpu time    0.7829: real time    0.7848 \n\
       DOS:  cpu time    0.0000: real time    0.0006 \n\
    -------------------------------------------- \n\
      LOOP:  cpu time    0.8239: real time    0.8272 \n\
 \n\
 eigenvalue-minimisations  :   768 \n\
 total energy-change (2. order) :-0.3365529E-05  (-0.1528204E-05) \n\
 number of electron      63.9999999 magnetization  \n\
 augmentation part       12.9708022 magnetization  \n\
 \n\
 Free energy of the ion-electron system (eV) \n\
  --------------------------------------------------- \n\
  alpha Z        PSCENC =       481.10416195 \n\
  Ewald energy   TEWEN  =     -5024.90579447 \n\
  -Hartree energ DENC   =      -973.95171042 \n\
  -exchange      EXHF   =         0.00000000 \n\
  -V(xc)+E(xc)   XCENC  =       294.02494192 \n\
  PAW double counting   =      5229.91755155    -5226.59116318 \n\
  entropy T*S    EENTRO =        -0.00000000 \n\
  eigenvalues    EBANDS =      -577.14209049 \n\
  atomic energy  EATOM  =      5753.77896827 \n\
  Solvation  Ediel_sol  =         0.00000000 \n\
  --------------------------------------------------- \n\
  free energy    TOTEN  =       -43.76513486 eV \n\
 \n\
  energy without entropy =      -43.76513486  energy(sigma->0) =      -43.76513486 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
 \n\
 average (electrostatic) potential at core \n\
  the test charge radii are     1.1137  0.7215 \n\
  (the norm of the test charge is              1.0000) \n\
       1 -40.3960       2 -39.6195       3 -40.1085       4 -40.5569       5 -68.0189 \n\
       6 -67.8638       7 -67.7498       8 -68.1938 \n\
  \n\
  \n\
  \n\
 E-fermi :   5.0877     XC(G=0): -12.1852     alpha+bet :-18.7323 \n\
 \n\
 \n\
 k-point     1 :       0.0000    0.0000    0.0000 \n\
  band No.  band energies     occupation  \n\
      1     -32.2756      2.00000 \n\
      2     -31.8709      2.00000 \n\
      3     -31.6231      2.00000 \n\
      4     -31.1938      2.00000 \n\
      5     -15.7848      2.00000 \n\
      6     -15.1631      2.00000 \n\
      7     -14.7749      2.00000 \n\
      8     -13.3421      2.00000 \n\
      9     -13.1843      2.00000 \n\
     10     -13.0670      2.00000 \n\
     11     -12.8561      2.00000 \n\
     12     -12.5247      2.00000 \n\
     13     -12.4469      2.00000 \n\
     14     -12.3227      2.00000 \n\
     15     -12.2190      2.00000 \n\
     16     -12.0708      2.00000 \n\
     17     -11.4036      2.00000 \n\
     18      -9.6060      2.00000 \n\
     19      -9.3683      2.00000 \n\
     20      -8.4415      2.00000 \n\
     21       1.1094      2.00000 \n\
     22       2.2563      2.00000 \n\
     23       2.3870      2.00000 \n\
     24       2.6258      2.00000 \n\
     25       2.7483      2.00000 \n\
     26       3.1985      2.00000 \n\
     27       3.3306      2.00000 \n\
     28       3.6188      2.00000 \n\
     29       3.8068      2.00000 \n\
     30       4.4177      2.00000 \n\
     31       4.6566      2.00000 \n\
     32       4.7232      2.00000 \n\
     33      10.0784      0.00000 \n\
     34      10.6253      0.00000 \n\
     35      10.9035      0.00000 \n\
     36      11.3465      0.00000 \n\
     37      12.0585      0.00000 \n\
     38      12.3918      0.00000 \n\
     39      13.2211      0.00000 \n\
     40      13.7434      0.00000 \n\
     41      13.9470      0.00000 \n\
     42      13.9746      0.00000 \n\
     43      14.2283      0.00000 \n\
     44      14.6140      0.00000 \n\
     45      14.9210      0.00000 \n\
     46      15.1047      0.00000 \n\
     47      15.2777      0.00000 \n\
     48      15.5748      0.00000 \n\
 \n\
 k-point     2 :       0.3333    0.0000    0.0000 \n\
  band No.  band energies     occupation  \n\
      1     -32.0605      2.00000 \n\
      2     -31.8419      2.00000 \n\
      3     -31.5593      2.00000 \n\
      4     -31.1155      2.00000 \n\
      5     -15.3511      2.00000 \n\
      6     -14.8563      2.00000 \n\
      7     -14.6405      2.00000 \n\
      8     -14.5093      2.00000 \n\
      9     -13.7230      2.00000 \n\
     10     -13.4697      2.00000 \n\
     11     -13.2501      2.00000 \n\
     12     -13.1189      2.00000 \n\
     13     -12.8441      2.00000 \n\
     14     -12.7005      2.00000 \n\
     15     -12.4956      2.00000 \n\
     16     -12.2416      2.00000 \n\
     17      -9.6699      2.00000 \n\
     18      -9.5902      2.00000 \n\
     19      -9.2994      2.00000 \n\
     20      -9.0410      2.00000 \n\
     21       2.1327      2.00000 \n\
     22       2.2200      2.00000 \n\
     23       2.3708      2.00000 \n\
     24       2.5690      2.00000 \n\
     25       2.9911      2.00000 \n\
     26       3.1158      2.00000 \n\
     27       3.4173      2.00000 \n\
     28       3.4921      2.00000 \n\
     29       3.6123      2.00000 \n\
     30       3.8621      2.00000 \n\
     31       3.9728      2.00000 \n\
     32       4.2345      2.00000 \n\
     33      10.3438      0.00000 \n\
     34      10.8943      0.00000 \n\
     35      11.1744      0.00000 \n\
     36      11.9265      0.00000 \n\
     37      12.4639      0.00000 \n\
     38      13.2153      0.00000 \n\
     39      13.4389      0.00000 \n\
     40      13.7759      0.00000 \n\
     41      13.9742      0.00000 \n\
     42      14.5840      0.00000 \n\
     43      14.7530      0.00000 \n\
     44      15.0864      0.00000 \n\
     45      15.3697      0.00000 \n\
     46      15.6245      0.00000 \n\
     47      16.0311      0.00000 \n\
     48      16.1297      0.00000 \n\
 \n\
 k-point     3 :       0.0000    0.5000   -0.0000 \n\
  band No.  band energies     occupation  \n\
      1     -32.1721      2.00000 \n\
      2     -31.9775      2.00000 \n\
      3     -31.6038      2.00000 \n\
      4     -31.2193      2.00000 \n\
      5     -15.8516      2.00000 \n\
      6     -14.7658      2.00000 \n\
      7     -14.4578      2.00000 \n\
      8     -14.0692      2.00000 \n\
      9     -13.3217      2.00000 \n\
     10     -13.1257      2.00000 \n\
     11     -12.7373      2.00000 \n\
     12     -12.5771      2.00000 \n\
     13     -12.5152      2.00000 \n\
     14     -12.4639      2.00000 \n\
     15     -12.2091      2.00000 \n\
     16     -12.0490      2.00000 \n\
     17     -11.1131      2.00000 \n\
     18      -9.5725      2.00000 \n\
     19      -9.2904      2.00000 \n\
     20      -8.2831      2.00000 \n\
     21       1.1861      2.00000 \n\
     22       2.2839      2.00000 \n\
     23       2.4386      2.00000 \n\
     24       2.6673      2.00000 \n\
     25       2.9387      2.00000 \n\
     26       3.1689      2.00000 \n\
     27       3.2957      2.00000 \n\
     28       3.6766      2.00000 \n\
     29       3.8449      2.00000 \n\
     30       3.9227      2.00000 \n\
     31       4.2228      2.00000 \n\
     32       4.3596      2.00000 \n\
     33       8.3334      0.00000 \n\
     34      11.5109      0.00000 \n\
     35      11.8676      0.00000 \n\
     36      12.0032      0.00000 \n\
     37      12.3128      0.00000 \n\
     38      12.6716      0.00000 \n\
     39      12.9633      0.00000 \n\
     40      13.7347      0.00000 \n\
     41      13.8867      0.00000 \n\
     42      14.1648      0.00000 \n\
     43      14.5698      0.00000 \n\
     44      14.7494      0.00000 \n\
     45      15.2218      0.00000 \n\
     46      15.4653      0.00000 \n\
     47      15.6884      0.00000 \n\
     48      15.9075      0.00000 \n\
 \n\
 k-point     4 :       0.3333    0.5000   -0.0000 \n\
  band No.  band energies     occupation  \n\
      1     -32.0408      2.00000 \n\
      2     -31.8530      2.00000 \n\
      3     -31.5578      2.00000 \n\
      4     -31.1278      2.00000 \n\
      5     -15.3579      2.00000 \n\
      6     -14.8944      2.00000 \n\
      7     -14.7318      2.00000 \n\
      8     -14.5741      2.00000 \n\
      9     -13.6351      2.00000 \n\
     10     -13.4520      2.00000 \n\
     11     -13.1231      2.00000 \n\
     12     -12.9525      2.00000 \n\
     13     -12.7805      2.00000 \n\
     14     -12.6947      2.00000 \n\
     15     -12.5045      2.00000 \n\
     16     -12.4099      2.00000 \n\
     17      -9.8449      2.00000 \n\
     18      -9.4549      2.00000 \n\
     19      -9.2680      2.00000 \n\
     20      -9.1123      2.00000 \n\
     21       2.1147      2.00000 \n\
     22       2.2546      2.00000 \n\
     23       2.3949      2.00000 \n\
     24       2.4533      2.00000 \n\
     25       2.9907      2.00000 \n\
     26       3.1737      2.00000 \n\
     27       3.3035      2.00000 \n\
     28       3.4271      2.00000 \n\
     29       3.6631      2.00000 \n\
     30       3.8557      2.00000 \n\
     31       4.1059      2.00000 \n\
     32       4.3112      2.00000 \n\
     33       9.9191      0.00000 \n\
     34      10.3869      0.00000 \n\
     35      11.5488      0.00000 \n\
     36      11.7864      0.00000 \n\
     37      12.2486      0.00000 \n\
     38      13.0856      0.00000 \n\
     39      13.6892      0.00000 \n\
     40      14.2541      0.00000 \n\
     41      14.3916      0.00000 \n\
     42      14.5982      0.00000 \n\
     43      14.8632      0.00000 \n\
     44      15.2554      0.00000 \n\
     45      15.3942      0.00000 \n\
     46      15.7093      0.00000 \n\
     47      15.7883      0.00000 \n\
     48      15.9821      0.00000 \n\
 \n\
 k-point     5 :       0.0000    0.0000    0.5000 \n\
  band No.  band energies     occupation  \n\
      1     -32.2006      2.00000 \n\
      2     -31.9249      2.00000 \n\
      3     -31.6469      2.00000 \n\
      4     -31.2008      2.00000 \n\
      5     -15.2684      2.00000 \n\
      6     -15.0068      2.00000 \n\
      7     -14.4899      2.00000 \n\
      8     -14.1273      2.00000 \n\
      9     -13.5245      2.00000 \n\
     10     -13.1046      2.00000 \n\
     11     -13.0587      2.00000 \n\
     12     -12.7684      2.00000 \n\
     13     -12.6278      2.00000 \n\
     14     -12.5283      2.00000 \n\
     15     -12.3014      2.00000 \n\
     16     -12.0972      2.00000 \n\
     17      -9.7398      2.00000 \n\
     18      -9.5644      2.00000 \n\
     19      -9.3494      2.00000 \n\
     20      -9.0086      2.00000 \n\
     21       1.8715      2.00000 \n\
     22       2.0497      2.00000 \n\
     23       2.3653      2.00000 \n\
     24       2.6815      2.00000 \n\
     25       2.7381      2.00000 \n\
     26       3.1423      2.00000 \n\
     27       3.4814      2.00000 \n\
     28       3.7022      2.00000 \n\
     29       3.8204      2.00000 \n\
     30       4.0397      2.00000 \n\
     31       4.2397      2.00000 \n\
     32       4.2821      2.00000 \n\
     33      10.3487      0.00000 \n\
     34      10.5063      0.00000 \n\
     35      11.2662      0.00000 \n\
     36      11.9862      0.00000 \n\
     37      12.5972      0.00000 \n\
     38      12.6195      0.00000 \n\
     39      13.0591      0.00000 \n\
     40      13.3205      0.00000 \n\
     41      13.6772      0.00000 \n\
     42      13.8568      0.00000 \n\
     43      13.9662      0.00000 \n\
     44      14.2597      0.00000 \n\
     45      14.8838      0.00000 \n\
     46      14.9189      0.00000 \n\
     47      15.0695      0.00000 \n\
     48      15.3656      0.00000 \n\
 \n\
 k-point     6 :       0.3333    0.0000    0.5000 \n\
  band No.  band energies     occupation  \n\
      1     -32.0404      2.00000 \n\
      2     -31.8405      2.00000 \n\
      3     -31.5788      2.00000 \n\
      4     -31.1195      2.00000 \n\
      5     -15.5511      2.00000 \n\
      6     -14.9888      2.00000 \n\
      7     -14.6548      2.00000 \n\
      8     -14.2705      2.00000 \n\
      9     -13.7854      2.00000 \n\
     10     -13.3198      2.00000 \n\
     11     -13.1502      2.00000 \n\
     12     -12.9531      2.00000 \n\
     13     -12.7962      2.00000 \n\
     14     -12.5990      2.00000 \n\
     15     -12.4978      2.00000 \n\
     16     -12.2688      2.00000 \n\
     17     -10.5982      2.00000 \n\
     18      -9.4830      2.00000 \n\
     19      -9.1738      2.00000 \n\
     20      -8.6400      2.00000 \n\
     21       1.5815      2.00000 \n\
     22       2.1988      2.00000 \n\
     23       2.5256      2.00000 \n\
     24       2.6930      2.00000 \n\
     25       2.9386      2.00000 \n\
     26       3.0720      2.00000 \n\
     27       3.3922      2.00000 \n\
     28       3.5438      2.00000 \n\
     29       3.7012      2.00000 \n\
     30       3.8972      2.00000 \n\
     31       4.0373      2.00000 \n\
     32       4.2534      2.00000 \n\
     33       9.1831      0.00000 \n\
     34      10.6489      0.00000 \n\
     35      11.2632      0.00000 \n\
     36      12.0643      0.00000 \n\
     37      12.8914      0.00000 \n\
     38      13.4991      0.00000 \n\
     39      14.0103      0.00000 \n\
     40      14.2874      0.00000 \n\
     41      14.5799      0.00000 \n\
     42      14.7044      0.00000 \n\
     43      14.9470      0.00000 \n\
     44      15.2664      0.00000 \n\
     45      15.5223      0.00000 \n\
     46      15.6975      0.00000 \n\
     47      15.9553      0.00000 \n\
     48      16.1192      0.00000 \n\
 \n\
 k-point     7 :       0.0000    0.5000    0.5000 \n\
  band No.  band energies     occupation  \n\
      1     -32.1696      2.00000 \n\
      2     -31.8995      2.00000 \n\
      3     -31.6694      2.00000 \n\
      4     -31.2433      2.00000 \n\
      5     -15.2008      2.00000 \n\
      6     -14.8676      2.00000 \n\
      7     -14.7313      2.00000 \n\
      8     -14.4658      2.00000 \n\
      9     -13.1473      2.00000 \n\
     10     -13.0493      2.00000 \n\
     11     -12.8642      2.00000 \n\
     12     -12.7508      2.00000 \n\
     13     -12.6003      2.00000 \n\
     14     -12.4804      2.00000 \n\
     15     -12.4701      2.00000 \n\
     16     -12.0413      2.00000 \n\
     17     -10.1084      2.00000 \n\
     18      -9.5760      2.00000 \n\
     19      -9.1654      2.00000 \n\
     20      -8.9923      2.00000 \n\
     21       1.8793      2.00000 \n\
     22       2.1240      2.00000 \n\
     23       2.3075      2.00000 \n\
     24       2.5858      2.00000 \n\
     25       2.6550      2.00000 \n\
     26       3.2843      2.00000 \n\
     27       3.3687      2.00000 \n\
     28       3.5120      2.00000 \n\
     29       3.8613      2.00000 \n\
     30       4.0571      2.00000 \n\
     31       4.2630      2.00000 \n\
     32       4.4605      2.00000 \n\
     33       9.2488      0.00000 \n\
     34       9.5165      0.00000 \n\
     35      11.0569      0.00000 \n\
     36      11.6670      0.00000 \n\
     37      12.6272      0.00000 \n\
     38      13.1384      0.00000 \n\
     39      13.6403      0.00000 \n\
     40      13.9521      0.00000 \n\
     41      14.0007      0.00000 \n\
     42      14.7023      0.00000 \n\
     43      14.8539      0.00000 \n\
     44      15.0156      0.00000 \n\
     45      15.0968      0.00000 \n\
     46      15.1509      0.00000 \n\
     47      15.5567      0.00000 \n\
     48      15.7374      0.00000 \n\
 \n\
 k-point     8 :       0.3333    0.5000    0.5000 \n\
  band No.  band energies     occupation  \n\
      1     -32.0291      2.00000 \n\
      2     -31.8427      2.00000 \n\
      3     -31.5753      2.00000 \n\
      4     -31.1338      2.00000 \n\
      5     -15.6140      2.00000 \n\
      6     -15.0717      2.00000 \n\
      7     -14.7282      2.00000 \n\
      8     -14.3957      2.00000 \n\
      9     -13.4236      2.00000 \n\
     10     -13.1806      2.00000 \n\
     11     -13.0345      2.00000 \n\
     12     -12.9506      2.00000 \n\
     13     -12.8147      2.00000 \n\
     14     -12.6469      2.00000 \n\
     15     -12.4735      2.00000 \n\
     16     -12.3394      2.00000 \n\
     17     -10.4440      2.00000 \n\
     18      -9.6426      2.00000 \n\
     19      -9.4265      2.00000 \n\
     20      -8.5805      2.00000 \n\
     21       1.5823      2.00000 \n\
     22       2.0391      2.00000 \n\
     23       2.4958      2.00000 \n\
     24       2.6590      2.00000 \n\
     25       2.9540      2.00000 \n\
     26       3.1288      2.00000 \n\
     27       3.3629      2.00000 \n\
     28       3.5190      2.00000 \n\
     29       3.7981      2.00000 \n\
     30       4.0937      2.00000 \n\
     31       4.2091      2.00000 \n\
     32       4.4096      2.00000 \n\
     33       8.6172      0.00000 \n\
     34       9.7763      0.00000 \n\
     35      11.6969      0.00000 \n\
     36      12.1906      0.00000 \n\
     37      12.6320      0.00000 \n\
     38      13.3254      0.00000 \n\
     39      14.2496      0.00000 \n\
     40      14.4613      0.00000 \n\
     41      14.8719      0.00000 \n\
     42      14.9598      0.00000 \n\
     43      15.2352      0.00000 \n\
     44      15.3797      0.00000 \n\
     45      15.4576      0.00000 \n\
     46      15.7815      0.00000 \n\
     47      15.8860      0.00000 \n\
     48      16.1355      0.00000 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 soft charge-density along one line, spin component           1 \n\
         0         1         2         3         4         5         6         7         8         9 \n\
 total charge-density along one line \n\
  \n\
 pseudopotential strength for first ion, spin component:           1 \n\
-17.703  -3.036  -0.032  -0.001   0.007   0.045   0.001  -0.010 \n\
 -3.036  -3.685  -0.018   0.000   0.005  -0.002  -0.000   0.000 \n\
 -0.032  -0.018 -17.652  -0.012   0.007   4.637   0.011  -0.006 \n\
 -0.001   0.000  -0.012 -17.639   0.008   0.011   4.624  -0.008 \n\
  0.007   0.005   0.007   0.008 -17.640  -0.006  -0.008   4.626 \n\
  0.045  -0.002   4.637   0.011  -0.006   9.264  -0.006   0.003 \n\
  0.001  -0.000   0.011   4.624  -0.008  -0.006   9.270   0.004 \n\
 -0.010   0.000  -0.006  -0.008   4.626   0.003   0.004   9.270 \n\
  0.006   0.001   0.008  -0.001  -0.030  -0.007  -0.001   0.029 \n\
 -0.011  -0.002  -0.000  -0.028  -0.001   0.002   0.030  -0.001 \n\
  0.004   0.001   0.016  -0.003  -0.003  -0.017  -0.001   0.001 \n\
  0.008   0.002  -0.001   0.007   0.002  -0.001  -0.010   0.001 \n\
  0.005   0.001   0.024   0.001   0.004  -0.029  -0.001  -0.005 \n\
  0.008   0.001   0.010  -0.001  -0.037  -0.008  -0.001   0.034 \n\
 -0.014  -0.002  -0.001  -0.034  -0.001   0.002   0.034  -0.001 \n\
  0.005   0.001   0.019  -0.004  -0.003  -0.019  -0.000   0.002 \n\
  0.009   0.002  -0.001   0.009   0.002  -0.001  -0.011   0.000 \n\
  0.007   0.001   0.030   0.001   0.005  -0.033  -0.001  -0.006 \n\
 total augmentation occupancy for first ion, spin component:           1 \n\
  2.007   0.036   0.002  -0.000  -0.001   0.002  -0.000  -0.001  -0.001   0.003  -0.001  -0.001  -0.001  -0.001   0.001  -0.001 \n\
  0.036   0.967   0.006  -0.003  -0.001   0.056   0.005  -0.021   0.035  -0.089   0.043   0.078   0.053  -0.021   0.055  -0.025 \n\
  0.002   0.006   2.015  -0.000   0.000   0.023   0.003  -0.002   0.005   0.015  -0.006  -0.012  -0.039   0.000  -0.006   0.002 \n\
 -0.000  -0.003  -0.000   2.015   0.000   0.003   0.020  -0.002  -0.015   0.019  -0.033  -0.026   0.002   0.005  -0.008   0.011 \n\
 -0.001  -0.001   0.000   0.000   2.016  -0.002  -0.002   0.020  -0.011  -0.012  -0.012   0.019  -0.009   0.003   0.003   0.003 \n\
  0.002   0.056   0.023   0.003  -0.002   0.068  -0.006   0.004  -0.051  -0.038  -0.024   0.037   0.041   0.031   0.023   0.015 \n\
 -0.000   0.005   0.003   0.020  -0.002  -0.006   0.074   0.004   0.042   0.031   0.103   0.037  -0.012  -0.025  -0.020  -0.062 \n\
 -0.001  -0.021  -0.002  -0.002   0.020   0.004   0.004   0.074   0.118   0.036   0.045  -0.061   0.010  -0.073  -0.022  -0.027 \n\
 -0.001   0.035   0.005  -0.015  -0.011  -0.051   0.042   0.118   1.162  -0.018  -0.008  -0.398  -0.052  -0.706   0.011   0.006 \n\
  0.003  -0.089   0.015   0.019  -0.012  -0.038   0.031   0.036  -0.018   0.721   0.384   0.018  -0.060   0.011  -0.430  -0.237 \n\
 -0.001   0.043  -0.006  -0.033  -0.012  -0.024   0.103   0.045  -0.008   0.384   1.005   0.130  -0.080   0.006  -0.238  -0.604 \n\
 -0.001   0.078  -0.012  -0.026   0.019   0.037   0.037  -0.061  -0.398   0.018   0.130   0.668   0.049   0.249  -0.012  -0.081 \n\
 -0.001   0.053  -0.039   0.002  -0.009   0.041  -0.012   0.010  -0.052  -0.060  -0.080   0.049   0.450   0.031   0.037   0.050 \n\
 -0.001  -0.021   0.000   0.005   0.003   0.031  -0.025  -0.073  -0.706   0.011   0.006   0.249   0.031   0.435  -0.007  -0.004 \n\
  0.001   0.055  -0.006  -0.008   0.003   0.023  -0.020  -0.022   0.011  -0.430  -0.238  -0.012   0.037  -0.007   0.259   0.149 \n\
 -0.001  -0.025   0.002   0.011   0.003   0.015  -0.062  -0.027   0.006  -0.237  -0.604  -0.081   0.050  -0.004   0.149   0.368 \n\
 -0.001  -0.046   0.004   0.010  -0.006  -0.022  -0.022   0.037   0.249  -0.011  -0.080  -0.398  -0.030  -0.157   0.007   0.051 \n\
 -0.001  -0.032   0.014  -0.000   0.003  -0.024   0.007  -0.005   0.031   0.037   0.050  -0.030  -0.261  -0.019  -0.023  -0.031 \n\
 \n\
 \n\
------------------------ aborting loop because EDIFF is reached ---------------------------------------- \n\
 \n\
 \n\
    CHARGE:  cpu time    0.0910: real time    0.0907 \n\
    FORLOC:  cpu time    0.0010: real time    0.0009 \n\
    FORNL :  cpu time    0.3080: real time    0.3096 \n\
    STRESS:  cpu time    0.8599: real time    0.8609 \n\
    FORCOR:  cpu time    0.0140: real time    0.0132 \n\
    FORHAR:  cpu time    0.0030: real time    0.0032 \n\
    MIXING:  cpu time    0.0010: real time    0.0009 \n\
    OFIELD:  cpu time    0.0000: real time    0.0001 \n\
 \n\
  FORCE on cell =-STRESS in cart. coord.  units (eV): \n\
  Direction    XX          YY          ZZ          XY          YZ          ZX \n\
  -------------------------------------------------------------------------------------- \n\
  Alpha Z   481.10416   481.10416   481.10416 \n\
  Ewald   -1688.43537 -1707.61433 -1628.86280    70.52319    -7.58321    -1.47527 \n\
  Hartree   318.50467   305.83100   349.60230    38.26781    -4.64140    -0.33097 \n\
  E(xc)    -306.43642  -306.50035  -306.24475     0.22523    -0.06017     0.03500 \n\
  Local     304.10936   337.51224   212.38500  -110.50078    13.57706     0.51777 \n\
  n-local  -275.20550  -274.70869  -275.79227    -0.78707    -0.05748    -0.11247 \n\
  augment   123.95945   123.57237   124.39682     0.70140    -0.26112     0.20038 \n\
  Kinetic  1082.92481  1081.03648  1085.36424     4.46528     0.08303     0.03524 \n\
  Fock        0.00000     0.00000     0.00000     0.00000     0.00000     0.00000 \n\
  ------------------------------------------------------------------------------------- \n\
  Total      40.52516    40.23287    41.95271     2.89506     1.05670    -1.13033 \n\
  in kB     812.88607   807.02309   841.52106    58.07145    21.19622   -22.67300 \n\
  external pressure =      170.48 kB  Pullay stress =      650.00 kB \n\
 \n\
 \n\
 VOLUME and BASIS-vectors are now : \n\
 ----------------------------------------------------------------------------- \n\
  energy-cutoff  :      500.00 \n\
  volume of cell :       79.87 \n\
      direct lattice vectors                 reciprocal lattice vectors \n\
     3.063970000  0.000000000  0.000000000     0.326373953 -0.088186988 -0.038529986 \n\
     1.355640000  5.017130000  0.000000000     0.000000000  0.199317139  0.043748454 \n\
     0.305250000 -1.140470000  5.195960000     0.000000000  0.000000000  0.192457217 \n\
 \n\
  length of vectors \n\
     3.063970000  5.197052361  5.328400295     0.340266751  0.204061876  0.192457217 \n\
 \n\
 \n\
 FORCES acting on ions \n\
    electron-ion (+dipol)            ewald-force                    non-local-force                 convergence-correction \n\
 ----------------------------------------------------------------------------------------------- \n\
   -.118E+02 0.831E+02 -.658E+01   0.133E+02 -.883E+02 0.676E+01   -.251E+00 0.968E+00 -.259E+00   -.992E-03 0.523E-02 -.351E-02 \n\
   0.116E+02 -.385E+02 0.453E+00   -.119E+02 0.409E+02 0.157E+00   -.369E+00 -.119E+01 -.438E+00   0.998E-03 -.244E-02 -.731E-03 \n\
   -.154E+02 0.199E+02 0.111E+02   0.138E+02 -.187E+02 -.124E+02   0.187E+01 -.124E+01 0.970E+00   -.554E-03 -.548E-04 0.296E-02 \n\
   -.169E+02 -.144E+01 -.155E+02   0.168E+02 0.365E+01 0.151E+02   0.682E+00 -.194E+01 0.508E+00   -.131E-02 0.174E-02 -.997E-04 \n\
   0.639E+01 -.678E+01 0.110E+02   -.626E+01 0.629E+01 -.102E+02   -.342E+00 0.146E+01 -.103E+01   -.473E-03 0.229E-02 0.233E-02 \n\
   0.405E+01 -.168E+02 -.214E+02   -.368E+01 0.165E+02 0.234E+02   -.136E+01 0.162E+01 -.181E+01   -.668E-03 -.376E-03 0.144E-02 \n\
   0.192E+02 -.354E+02 0.166E+02   -.209E+02 0.377E+02 -.185E+02   0.148E+01 -.192E+01 0.191E+01   -.108E-02 -.826E-04 -.183E-02 \n\
   0.450E+00 -.545E-01 0.284E+01   -.106E+01 0.193E+01 -.437E+01   0.737E+00 -.188E+01 0.162E+01   0.541E-04 0.263E-02 -.240E-02 \n\
 ----------------------------------------------------------------------------------------------- \n\
   -.245E+01 0.411E+01 -.147E+01   -.371E-13 0.320E-13 0.711E-13   0.245E+01 -.412E+01 0.147E+01   -.402E-02 0.894E-02 -.184E-02 \n\
  \n\
  \n\
 POSITION                                       TOTAL-FORCE (eV/Angst) \n\
 ----------------------------------------------------------------------------------- \n\
      1.14652      3.32842      1.59978         1.168281     -4.245817     -0.080442 \n\
      2.96550      0.81050      0.27757        -0.735246      1.188434      0.171096 \n\
      3.10346      2.41291      4.21309         0.238719      0.022373     -0.258197 \n\
      1.42766      0.32431      2.90111         0.595000      0.266135      0.121868 \n\
      1.55942      3.72053      3.58402        -0.209437      0.968435     -0.291308 \n\
      1.75007      0.96799      4.89189        -0.990043      1.383871      0.212539 \n\
      2.95605      1.63997      2.22538        -0.192415      0.419208      0.024764 \n\
      2.72697      4.51168      0.91693         0.125140     -0.002640      0.099680 \n\
 ----------------------------------------------------------------------------------- \n\
    total drift:                                0.000484     -0.000437      0.000143 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 \n\
  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV) \n\
  --------------------------------------------------- \n\
  free  energy   TOTEN  =       -43.76513486 eV \n\
 \n\
  energy  without entropy=      -43.76513486  energy(sigma->0) =      -43.76513486 \n\
  enthalpy is  TOTEN    =       -11.36040269 eV   P V=       32.40473217 \n\
 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
    POTLOK:  cpu time    0.0410: real time    0.0417 \n\
 \n\
 \n\
-------------------------------------------------------------------------------------------------------- \n\
 \n\
 \n\
 writing wavefunctions \n\
     LOOP+:  cpu time   13.0910: real time   13.1306 \n\
    4ORBIT:  cpu time    0.0000: real time    0.0000 \n\
 \n\
 total amount of memory used by VASP MPI-rank0    38421. kBytes \n\
======================================================================= \n\
 \n\
   base      :      30000. kBytes \n\
   nonl-proj :       4980. kBytes \n\
   fftplans  :       1298. kBytes \n\
   grid      :       1190. kBytes \n\
   one-center:        124. kBytes \n\
   wavefun   :        829. kBytes \n\
  \n\
   \n\
   \n\
 General timing and accounting informations for this job: \n\
 ======================================================== \n\
   \n\
                  Total CPU time used (sec):       15.618 \n\
                            User time (sec):       15.082 \n\
                          System time (sec):        0.536 \n\
                         Elapsed time (sec):       15.578 \n\
   \n\
                   Maximum memory used (kb):       80212. \n\
                   Average memory used (kb):           0. \n\
   \n\
                          Minor page faults:        13072 \n\
                          Major page faults:            0 \n\
                 Voluntary context switches:          331 \n\
";
