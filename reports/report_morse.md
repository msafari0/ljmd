# Report on implementation of Morse potential as an interatomic potential in the system of Argon.
date: 03/04/2023 - 07/04/2023
student: Mandana Safari

## Introduction:
The purpose of this report is to discuss the implementation of Morse potential as an interatomic potential in the system of Argon.
The goal of the project was to understand the behavior of Argon atoms and molecules in different forces by simulating their interactions (here using Morse potential).
In this report, I will discuss the implementation process and the results obtained from the simulations.

## Background:
Morse potential is an interatomic potential that is used to describe the interactions between atoms and molecules as Van del Waals interactions. It is a more accurate model 
than other potentials like the lenard jones potential, as it takes into account the anharmonicity of the interatomic interactions. Therefore, it is more used in various
scientific fields like physics, chemistry, and materials science.

## Methodology:
 The implementation process involved considering a shift in potential energy in a way that energy and force at cutoff becomes zero. There are three parameters for Morse
 potential energy, which can be derived from lenard jones constant as well.

To implement Morse potential, I am adding to the existing Molecular dynamic code to include the Morse potential equations for calculating the energy and force between atoms
by considering a force flag `fflag` in input file `argon_108.inp`. I considered 0 for lenard jones potential and 1 for Morse potential. In this way, we can compile once and
based on needs use different force fields.
We also added a shift in potential energy to ensure that the energy and force at cutoff became zero.
Morse potential is called by this formula: `V(r) = D*(1-exp(-a*(r-r0)))^2`. Although this version is used for investigation of bond features, here I used the smooth version:
`V(r) = D * (exp(-a*dr)*exp(-a*dr) - 2.0*exp(-a*dr)) - phi_rcut`. This formula comes from taylor expantion of `phi = D * (exp(-a*dr)*exp(-a*dr) - 2.0*exp(-a*dr))` around r_cut.
(Explanation based on https://docs.lammps.org/pair_morse.html)



```void morse_force(mdsys_t *sys)
{

    double a, r0,rcut, D, rcsq, rsq;
    double dr, r,rx,ry,rz;
    int i,j;

    ...
    
    rcsq=sys->rcut*sys->rcut; //square of the cutoff radius
    D = sys->epsilon;
    r0 = sys->sigma*pow(2.0,1.0/6.0);
    a = 2.0;

    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) { 

            ...
            
            double rsq = rx*rx + ry*ry + rz*rz;  //remove the expensive sqrt() function from the inner loop

            /* compute force and energy if within cutoff */
           if (rsq < rcsq) {
                /*force between 2 atoms */   
                r = sqrt(rsq); //remove the expensive sqrt() function from the inner loop
                dr = r - r0;   // define dr
                double dexp = exp(-a*dr); 
                double morse1 = 2.0*D*a; //introduce morse1 coefficient
                double alpha_dr = -a * (rcut - r0);  //introduce alpha_dr
                double phi_rcut  = D * (exp(2.0*alpha_dr) - 2.0*exp(alpha_dr)); //the energy at the cutoff radius
                double der_at_cutoff = -2.0*a*D * (exp(2.0*alpha_dr) - exp(alpha_dr)); //the derivative at the cutoff radius

                double fpartial = morse1*(dexp*dexp - dexp) / r;  //introduce fpartial, which is the force at the cutoff radius
                double ffac = D * ( fpartial + der_at_cutoff / r);  
        
                /*morse energy pair wise*/  
                sys->epot = D * (dexp*dexp - 2.0*dexp) - phi_rcut;  // the energy at the cutoff radius is subtracted to avoid double counting
                dr = rcut - r;  // redefine dr, since the energy at the cutoff radius is subtracted
                sys->epot += dr * der_at_cutoff; //adding the energy at the cutoff radius to the energy at the current distance
                
                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                sys->fx[j] -= rx*ffac;  
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;
            }
        }
    }
} ```

The gtest fucntion was coded accordingly.


## Results:
The Morse potential was a more accurate model than lenard jones potential, as it took into account the anharmonicity of the interatomic interactions.
It needs some consideration based of taylor expantion and a shift in energy and force. Without the shift, the simulations showed that the system was unstable,
and the atoms and molecules were not able to maintain their positions. However, with the shift, the system was stable, and the atoms and molecules were able to 
maintain their positions.

## Conclusion:
The simulations showed that Morse potential is a more accurate model than lenard jones potential, as it takes into account the anharmonicity of 
the interatomic interactions.The shift in potential energy is also essential for ensuring the stability of the system. This project has provided 
us with valuable insights into the behavior of atoms and molecules, which can be used to develop more accurate models for simulating their interactions.
