# Development

This package, **nap** , is an open source project. And you can modify
and use for free of charge and by your own risk.

Usually users have to modify program source codes when they apply the MD
program to specific subjects they are tackling. The source codes in this
package contains a lot of lines, but the algorithms usesd in this
package are mostly basic and users can easily understand most of these
if they are familiar with the molecular dynamics.

*Enjoy coding and your MD simulation ;)*

------------------------------------------------------------------------

## Formulas

For those who want to read and modify program sources, it would be
helpful to show some details and formulas of what the **pmd** is doing
in the program. Basically the **pmd** performs standard MD with the
velocity Verlet algorithm, one should firstly understand the standard MD
formulas using some standard textbooks[^Frenkel].

However, of course, there are a lot code-specific definitions to
implement the standard MD. This section shows some of the code-specific
definitions used in the **pmd**.

------------------------------------------------------------------------

### Simulation cell

The simulation cell consists of one scalar value, $l$, and three
vectors, $\mathbf{a}, \mathbf{b}, \mathbf{c}$, and is shown in the
atom-configuration file as, :

    5.472                     <--- l
    0.500    0.500   0.000    <--- ax, ay, az
    0.500    0.000   0.500    <--- bx, by, bz
    0.000    0.500   0.500    <--- cx, cy, cz

The cell information actually used in the code is the matrix,
$\mathbf{h}$, defined as,

$$\begin{aligned}
\begin{equation}
   \mathbf{h} = l \times (\mathbf{a},\mathbf{b},\mathbf{c})= l \times
   \begin{pmatrix}
   a_x & b_x & c_x \\
   a_y & b_y & c_y \\
   a_z & b_z & c_z
   \end{pmatrix}.
\end{equation}
\end{aligned}$$

------------------------------------------------------------------------

### Positions

Atom positions, which are defined as `ra(1:3,1:natm)` in the code, are
normalized within 0.0 and 1.0 during the MD simulation so that they become
absolute positions after multiplying the cell matrix, $\mathbf{h}$, as,

$$\begin{aligned}
\begin{eqnarray*}
\mathbf{r}' & = & \mathbf{h} \cdot \mathbf{r}, \\
            & = & r_a\mathbf{a} +r_b\mathbf{b} +r_c\mathbf{c},
\end{eqnarray*}
\end{aligned}$$

which is written in the Fortran code as follows,

```fortran
xi(1:3)= h(1:3,1)*ra(1,i) +h(1:3,2)*ra(2,i) +h(1:3,3)*ra(3,i)
```

------------------------------------------------------------------------

### Velocities

Atom velocities, `va(1:3,1:natm)`, are also scaled by the cell matrix.
However, not only that, but also scale by time, which means the
velocities in the code have actually length scale but velocity scale.
Thus, in the velocity Verlet algorithm, velocities are directly added to
the positions without multiplying time interval, $\Delta t$.

------------------------------------------------------------------------

### Accelerations

Atom accelerations, `aa(1:3,1:natm)`, are also scaled by the cell matrix
and are multplied by $\Delta t^2$ so that they become the units of
positions at the end of the force calculation subroutines. Thus, the
accelerations as well can be directly added to the positions without
multiplying $\Delta t^2$ in the velocity Verlet loop.

------------------------------------------------------------------------

### Velocity Verlet {#vv}

Basically MD with the velocity Verlet algorithm is very simple:

1.  Compute initial forces.
2.  Velocity Verlet loop starts.
    1.  Update velocities using the current forces with a half of time
        interval, $\Delta t/2$.
    2.  Update positions using the current velocities with a time
        interval, $\Delta t$.
    3.  Compute forces from the current positions.
    4.  Update velocities using the current forces with a half of time
        interval, $\Delta t/2$.

This is written mathematically as,

$$\begin{aligned}
\begin{eqnarray*}
\mathbf{v}_i^{(n+1)} &=& \mathbf{v}_i^* +\frac{\mathbf{f}_i^{(n)}}{m_i} \frac{\Delta t}{2}, \\
\mathbf{r}_i^{(n+1)} &=& \mathbf{r}_i^{(n)} +\mathbf{v}_i^{(n+1)}\Delta t, \\
\text{Compute}\ &\ & \mathbf{f}_i^{(n+1)}\left(\left\{\mathbf{r}^{(n+1)}\right\}\right), \\
\mathbf{v}_i^* &=& \mathbf{v}_i^{(n+1)} +\frac{\mathbf{f}_i^{(n+1)}}{m_i} \frac{\Delta t}{2},
\end{eqnarray*}
\end{aligned}$$

where subscript *i* is an atomic index and superscript *n* is an
MD-step.

This is implemented in the code as,

```Fortran
call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
&     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
&     ,mpi_md_world,myid_md,epi,epot0,nismax,acon,avol
&     ,cforce)
...
do istp=1,nstp
   ...
   va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)
   ...
   ra(1:3,1:natm)=ra(1:3,1:natm) +va(1:3,1:natm)
   ...
   call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
&         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
&         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
&         ,cforce)
   ...
   va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)
   ...
enddo
```

As you can see, the basic construction is very simple. However, the
codes for miscellaneous stuff such as thermostat, isobaric, and
parallelization hidden in the above code-block as `...` are rather
lengthy.

------------------------------------------------------------------------

### Berendsen barostat

Berendsen barostat is similar to Berendsen thermostat, which control the
stress or temperature moderately to the target ones.

------------------------------------------------------------------------

## Force implementation


Implementation of the interatomic force and potential is the core of MD
programming. And if you want to perform simulation that includes a
combination of elements that is not implemented in the **pmd**, you have
to implement the force calculation routine by yourself. Every force
routine is defined in a separated module such as `force_SW_Si.F90`,
which indicates the force routine of Stillinger-Weber type for Si
system, and is called via `get_force` subroutine in `force_common.F`
file. Thus you can write your own force routine by following `get_force`
and `force_SW_Si` subroutines.

The **force-type** used for the simulation is determined by a variable,
`cforce`, which is a character and specified in *in.pmd* file. The
corresponding force routine is called according to the **force-type**
name. Each force subroutine has to have the same arguments and order and
should be provided using `use` in the beginning of the `get_force`
routine.

------------------------------------------------------------------------

## ASE interface

There is a python script that connects pmd to **ASE (atomistic
simulation environment)**. This enable us small calculation of pmd much
easier. The following code shows how to use
`nappy/interface/ase/pmdrun.py` with ase program.

``` {.python}
import os,sys
from ase.io import read
from nappy.interface.ase.pmdrun import PMD

atoms=read('POSCAR',index=0,format='vasp')
os.system('cp /path/to/in.*.NN ./')
calc= PMD(label='pmd',command='/path/to/nap/pmd/pmd > out.pmd',force_type='NN')
atoms.set_calculator(calc)
print atoms.get_potential_energy()
print atoms.get_forces()
```

When `atoms.get_potential_energy()` is called, pmd program is performed
background and ASE gets calculation results from `out.pmd` file.

------------------------------------------------------------------------

## Versioning and tagging

The **nap** uses Git for version controlling. Developers are recommended
to make new branches to modify codes, and merge the change to the master
branch. Tags are named like `rev170605` which means the *revision* of
the date 2017-06-05. If you need to make new tag at the same date, which
would not happen often, you can name the tags like `rev170605_1` or
something like this by adding some suffix after an underscore.

Versioning in **nap** follows (not strictly) the [semantic
versioning](https://semver.org). However, since API is not defined, the
major version will never get one and the version will be like `v0.x,x`.

[^Frenkel]: Dann Frenkel and Berend Smit, *Understanding Molecular Simulation*, Academic Press
