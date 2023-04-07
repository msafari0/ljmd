## Report on individual task of implementing force optimization and reduction of expensive operators
date: 03/04/2023 - 07/04/2023
student: Mandana Safari

## Introduction:
The purpose of this report is to discuss the implementation of force optimization and the reduction of expensive operators in optimizing 
and parallelizing Molecular dynamic code with lenard jones potential. 
The goal of the project was to improve the performance of the code by optimizing it for parallel processing and reducing the computational
time required for each simulation. In this report, I will discuss the individual task of implementing force optimization and reducing the 
computational time by implementing the third law of newton.

## Methodology:
In this project, our group worked on optimizing and parallelizing Molecular dynamic code governed by lenard jones potential.
 My individual task was to optimize the force computation: including refactorization the code for better optimization and to avoid costly
operations or redundant work. Adapt data structures as needed. 
First, I implementted the third law of newton:

```static void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;

    ...
    
    for(i=0; i < (sys->natoms); ++i) {    
        for(j=0; j < (sys->natoms); ++j) {  ----------------------------->  "for(j=i+1; j < (sys->natoms); ++j) { "

            if (i==j) continue;

            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box); 
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            r = sqrt(rx*rx + ry*ry + rz*rz);

            if (r < sys->rcut) {
                ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);

                sys->epot += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)  
                                               -pow(sys->sigma/r,6.0));     ----> Omit the 0.5 coefficient

                sys->fx[i] += rx/r*ffac;   |                  | sys->fx[i] -= rx/r*ffac; sys->fy[i] -= ry/r*ffac; sys->fz[i] -= rz/r*ffac;
                sys->fy[i] += ry/r*ffac;   |----------------> | sys->fx[j] -= rx/r*ffac; sys->fy[j] -= ry/r*ffac; sys->fz[j] -= rz/r*ffac;
                sys->fz[i] += rz/r*ffac;   |
            }
        }
    }
}```

The original code was j=0, which means that the force on particle i was computed twice, then we needed a 0.5 coefficient in potential part.
This is fixed by starting the loop at j=i+1, and removing the coefficient.  The third law of newton states that for every action, there is
an equal and opposite reaction. Therefore, when two atoms interact, the force exerted on one atom is equal in magnitude but opposite in 
direction to the force exerted on the other atom. This is implementted with considering negative force components of J atoms at the end.

I also reduced the computational time by avoiding expensive operators like power and square roots. Instead, I used precomputed variables
to calculate these values, resulting in a further reduction in computational time.

```void force(mdsys_t *sys)
{
    double r,ffac;
    double rx,ry,rz;
    int i,j;
    //new variables for the optimization
    double c12,c6,rcsq, rsq, rm6, rm2;
    ...

    //remove the expensive sqrt(), pow() and division functions from the inner loop
    c12=4.0*sys->epsilon*pow(sys->sigma,12.0); //c12 is the 12th power of sigma
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0); //c6 is the 6th power of sigma
    rcsq=sys->rcut*sys->rcut; //square of the cutoff radius


    for(i=0; i < (sys->natoms); ++i) {
        for(j=i+1; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
           // if (i==j) continue; //it will be useless, since j starts from i+1

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz;  //remove the expensive sqrt() function from the inner loop

            /* compute force and energy if within cutoff */
           if (rsq < rcsq) {
                //remove the expensive pow() and division functions from the inner loop
                double rm6,rm2; rm2=1.0/rsq; rm6=rm2*rm2*rm2;
                ffac = (12.0*c12*rm6 - 6.0*c6)*rm6*rm2;
                sys->epot += rm6*(c12*rm6 - c6); 

                sys->fx[i] += rx*ffac;
                sys->fy[i] += ry*ffac;
                sys->fz[i] += rz*ffac;
                sys->fx[j] -= rx*ffac;  
                sys->fy[j] -= ry*ffac;
                sys->fz[j] -= rz*ffac;
            }
        }
    }
}```

Which is mostly by adding variables inwhich costly operations are precomputed once before loop. 

## Results:
After implementing force optimization and reducing the computational time by implementing the third law of newton and avoiding 
costly operations, I was able to achieve about 3 times speedup in total. 

                   | code type      | size | time     | speedup |
                   | ---------------|------|----------|---------|
                   | vanilla        | 108  | 1.562s   |  ....   |
                   | only force     | 108  | 0.793s   |  1.97   |
                   | force + OP opt | 108  | 0.523s   |  2.98   |
                   | vanilla        | 2916 | 294.744s |  ....   |
                   | only force     | 2916 | 143.956s |  2.04   |
                   | force + OP opt | 2916 | 147.796s |  1.99   |
                   
As we can see, with scaling the size of the system, effects of force and operations optimization are reduced, due to increasing 
the number of calculation between pairs.

## Conclusion:
In conclusion, implementing force optimization and reducing the computational time by implementing the third law of newton and 
avoiding expensive operators have improved the performance of the Molecular dynamic code about 3 times. 
