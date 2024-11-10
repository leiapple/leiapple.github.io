---
layout: post
title: Interstitial energy in LAMMPS
categories: Atomistic
description: interstitial energy
keywords: lammps, atomistic
---

# Calculating Interstitial Energy of Carbon in Iron Using LAMMPS

**Interstitial energy** refers to the energy required to insert an atom into the interstitial site of a crystal lattice. In this guide, we will calculate the interstitial energy of carbon in a body-centered cubic (BCC) iron lattice using the LAMMPS molecular dynamics software.

We will:
1. Set up a perfect BCC iron lattice and calculate its potential energy.
2. Insert a carbon atom into an octahedral interstitial site, relax the system and calculate its potential energy.
3. Calculate the potential energy of a carbon atom in diamond structure.
4. Compute the interstitial energy of carbon.

---

## Step 1: Setting Up the Simulation

First, install LAMMPS if you haven't already. You can download it from the [official LAMMPS website](https://www.lammps.org).

We will use an embedded-atom method (EAM) potential for the iron and carbon interaction. Make sure you have the potential files in your working directory. You can download the Fe-C EAM potential file from [OpenKIM](https://openkim.org) or other sources.

### Lattice Setup

We begin by defining a perfect BCC iron lattice. Below is the LAMMPS input script to set up the simulation box and create atoms for the BCC lattice. Here is the lammps command:

```
# Interstitial Energy of Carbon in BCC Iron

units metal
dimension 3
boundary p p p
atom_style atomic

# Define the BCC lattice with lattice constant 2.87 Angstroms (typical for Fe, but changes according to the potential)
lattice bcc 2.87
region box block 0 10 0 10 0 10
create_box 2 box

# Create Fe atoms in the BCC structure
create_atoms 1 box

# Define the potential (Fe-C alloy EAM potential)
pair_style eam
pair_coeff * * Fe_C.eam.alloy Fe C

# Assign masses to the atom types
mass 1 55.845  # Iron
mass 2 12.011  # Carbon

# Set up neighbor lists and thermodynamics output
neighbor 2.0 bin
neigh_modify delay 5

thermo 100
thermo_style custom step temp pe

# Minimize the system to get the energy of pure BCC iron
minimize 1.0e-6 1.0e-8 1000 10000

# Write the minimized configuration to a file
write_data bcc_iron.data

```

## Step 2: Insert Carbon Atom into Octahedral Interstitial Site

Once you’ve generated the BCC iron lattice and minimized the system, the next step is to insert a carbon atom into an octahedral interstitial site.

The octahedral interstitial site in a BCC unit cell is located at fractional coordinates such as (0.5, 0.5, 0), (0.5, 0, 0.5), and (0, 0.5, 0.5). We will insert a carbon atom at (0.5, 0.5, 0), which corresponds to one of the octahedral sites in the BCC lattice.

Add the following lines to the LAMMPS script:

```
# Insert a carbon atom at an octahedral site (0.5, 0.5, 0)
create_atoms 2 single 1.435 1.435 0.0 unit box # Position in Angstroms based on lattice constant 2.87 Å

# Minimize the system again after inserting the carbon atom
minimize 1.0e-6 1.0e-8 1000 10000

# Write the minimized configuration with the carbon interstitial
write_data fe_c_interstitial.data
```

To extract these energies, you can modify the LAMMPS script to output the potential energies at each stage:
```
# Output the total energy of the system
variable E_Fe equal pe  # Energy of pure BCC iron
print "Energy of pure BCC iron: ${E_Fe} eV"

# After inserting the carbon atom
variable E_Fe_C equal pe  # Energy of BCC iron with carbon interstitial
print "Energy of BCC iron with carbon interstitial: ${E_Fe_C} eV"
```

## Step 3: Calculate the potential energy of a carbon atom in diamond structure

The carbon energy can be obtained from a separate calculation in LAMMPS of carbon atomic energy in diamond structure. 

```
# -----------------------------------------------------------
# LAMMPS script to find the equilibrium lattice constant
# of carbon in a diamond structure and save results to file
# -----------------------------------------------------------

log             eos_diamond.log           
label           loop_start

variable        a equal 3.5                      # Initial guess for the lattice constant
variable        da equal 0.005                    # Increment for the lattice constant
variable        steps equal 20                   # Number of steps for varying the lattice constant

variable        i loop 1 ${steps}
variable        latparam equal ${a}+${da}*${i}

clear
units           metal
dimension       3
boundary        p p p
atom_style      atomic
atom_modify     map yes

#-------Define geometry  (2d X 2d) ---------------------
lattice         diamond ${latparam}
region          box block 0 1 0 1 0 1 units lattice
create_box      2 box
create_atoms    2 box

#-------Define interatomic potential--------------------
pair_style      eam
pair_coeff      * * Fe_C.eam.alloy Fe C
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

#-------Define atom mass (no needed for moelcular statics) 
mass            1 55.845
mass            2 12.011

#-------Compute------------------------------------------
compute 	eng all pe/atom
compute 	new all temp
compute 	csym all centro/atom bcc
compute 	poten all pe
compute 	stress all stress/atom NULL

#------Relaxation----------------------------------------
thermo          100
thermo_style    custom step pe lx ly lz pxx pyy pzz pxy pxz pyz press

min_style       cg
minimize        1e-30 10e-12 100000000 1000000000

#-----Calculate the lattice parameter--------------------
variable        tmp equal "lx"
variable        LX equal ${tmp}
variable        tmp equal "ly"
variable        LY equal ${tmp}
variable        tmp equal "lz"
variable        LZ equal ${tmp}
variable        tmp equal "atoms"
variable        N equal ${tmp}
variable        V equal (${LX}*${LY}*${LZ})
variable        alat equal (${LX}+${LY}+${LZ})/3
variable        tmp equal "pe"
variable        pe0 equal ${tmp}
variable        teng equal ${pe0}/${N}
variable        vol  equal ${V}/${N}

#-----output to datafolder-------------------------------- 
print           "${vol} ${teng} ${alat}" append ./volume_diamond.dat

#-----jump to next loop-----------------------------------
next            i
jump            in.eos_C_dia loop_start
```

## Step 4: Compute the interstitial energy of carbon.

Calculating Interstitial Energy
Now that we have the system both with and without the carbon atom, we can calculate the interstitial energy.
The interstitial energy is given by:
$E_{inter} = E_{Fe+C} - E_{Fe} - E_{C}$
Where:
- $E_{Fe+C}$ is the energy of the iron lattice with the carbon atom.
- $E_{Fe}$ is the energy of the pure iron lattice.
- $E_{C}$ is the energy of an isolated carbon atom.

By following this process, you can calculate the interstitial energy of carbon in iron using LAMMPS. This value is useful in understanding how carbon atoms interact with the iron lattice, which is important for studying phenomena like carbon diffusion in steels.
Feel free to modify the script to use different lattice sizes or to insert multiple interstitial atoms. The same approach can be extended to calculate interstitial energies for other elements in different metals.

Happy simulating!

