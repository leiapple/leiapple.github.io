# Calculating Interstitial Energy of Carbon in Iron Using LAMMPS

**Interstitial energy** refers to the energy required to insert an atom into the interstitial site of a crystal lattice. In this guide, we will calculate the interstitial energy of carbon in a body-centered cubic (BCC) iron lattice using the LAMMPS molecular dynamics software.

We will:
1. Set up a perfect BCC iron lattice.
2. Insert a carbon atom into an octahedral interstitial site.
3. Calculate the energy of the system with and without the carbon atom.
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

# Define the BCC lattice with lattice constant 2.87 Angstroms (typical for Fe)
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
create_atoms 2 single 1.435 1.435 0.0  # Position in Angstroms based on lattice constant 2.87 Å

# Minimize the system again after inserting the carbon atom
minimize 1.0e-6 1.0e-8 1000 10000

# Write the minimized configuration with the carbon interstitial
write_data fe_c_interstitial.data
```

## Step 3: Calculating Interstitial Energy

Now that we have the system both with and without the carbon atom, we can calculate the interstitial energy.

The interstitial energy is given by:

\[
E_{\text{interstitial}} = E_{\text{Fe+C}} - E_{\text{Fe}} - E_{\text{C, isolated}}
\]

Where:
- \( E_{\text{Fe+C}} \) is the energy of the iron lattice with the carbon atom.
- \( E_{\text{Fe}} \) is the energy of the pure iron lattice.
- \( E_{\text{C, isolated}} \) is the energy of an isolated carbon atom.

To extract these energies, you can modify the LAMMPS script to output the potential energies at each stage:

```
# Output the total energy of the system
variable E_Fe equal pe  # Energy of pure BCC iron
print "Energy of pure BCC iron: ${E_Fe} eV"

# After inserting the carbon atom
variable E_Fe_C equal pe  # Energy of BCC iron with carbon interstitial
print "Energy of BCC iron with carbon interstitial: ${E_Fe_C} eV"

# Run a separate simulation for a single carbon atom (to get its isolated energy)
# Define a single carbon atom in a box
```
The carbon energy can be obtained from a separate calculation in LAMMPS of carbon atomic energy in diamond structure.


## Step 4: Post-Processing and Calculation


By following this process, you can calculate the interstitial energy of carbon in iron using LAMMPS. This value is useful in understanding how carbon atoms interact with the iron lattice, which is important for studying phenomena like carbon diffusion in steels.

Feel free to modify the script to use different lattice sizes or to insert multiple interstitial atoms. The same approach can be extended to calculate interstitial energies for other elements in different metals.

Happy simulating!

