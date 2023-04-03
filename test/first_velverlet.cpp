#include <gtest/gtest.h>
extern "C" { // since our lib is C code, but the test here is C++
#include <mdlib.h>
}

TEST(FirstVelVerlet, ComputeHalfTimeStep) {

	int i;
    mdsys_t sys;

    /* fill in simple numbers */
    sys.natoms  = 2;
    sys.mass    = 2.2;
    sys.box     = 10;
    sys.dt      = 2;

    /* allocate memory for the system */
    sys.rx=new double[2]; 	sys.ry=new double[2]; 	sys.rz=new double[2];
	sys.vx=new double[2]; 	sys.vy=new double[2]; 	sys.vz=new double[2];
	sys.fx=new double[2]; 	sys.fy=new double[2]; 	sys.fz=new double[2];

	/*initialize positions with a simple cubic lattice */
    for (i=0; i<sys.natoms; ++i) {
        sys.rx[i] = 1 + (i%2)*sys.box/4.0;			
        sys.ry[i] = 1 + ((i/2)%2)*sys.box/4.0;		
        sys.rz[i] = 1 + (i/4)*sys.box/4.0;			
    }
    /* initialize velocities to fixed nonzero values */
    for (i=0; i<sys.natoms; ++i) {
        sys.vx[i] = 3*((i+1)*1.0);					
        sys.vy[i] = 3*(2.0/(i-0.5));				
        sys.vz[i] = 3*(-1.0*(i+1)*(i+1));		
    }
    /*initialize forces to fixed nonzero values */
    for (i=0; i<sys.natoms; ++i) {
        sys.fx[i] = 1500*(i+1 -(i%2)*2); 			
        sys.fy[i] = 1500*(2*i+1 -(i%3)*4);			
        sys.fz[i] = 1500*(-1.0*(i+1)*(i+1));		
    }

	/* propagate system with half step velocity verlet */
    velverlet_first_half(&sys);

	/* check that positions have been updated correctly */
	EXPECT_FLOAT_EQ(sys.rx[0], 7.570545) 	<< "rx[0] is not correct";
	EXPECT_FLOAT_EQ(sys.ry[0], -22.429455) 	<< "ry[0] is not correct";
	EXPECT_FLOAT_EQ(sys.rz[0], -5.570545)	<< "rz[0] is not correct";
	EXPECT_FLOAT_EQ(sys.rx[1], 15.500000)  	<< "rx[1] is not correct";
	EXPECT_FLOAT_EQ(sys.ry[1], 24.429455)  	<< "ry[1] is not correct";
	EXPECT_FLOAT_EQ(sys.rz[1], -25.282182) 	<< "rz[1] is not correct";

	/* check that velocities have been updated correctly */
	EXPECT_FLOAT_EQ(sys.vx[0], 3.285273) 	<< "vx[0] is not correct";
	EXPECT_FLOAT_EQ(sys.vy[0], -11.714727) 	<< "vy[0] is not correct";
	EXPECT_FLOAT_EQ(sys.vz[0], -3.285273)  	<< "vz[0] is not correct";
	EXPECT_FLOAT_EQ(sys.vx[1], 6.000000) 	<< "vx[1] is not correct";
	EXPECT_FLOAT_EQ(sys.vy[1], 11.714727) 	<< "vy[1] is not correct";
	EXPECT_FLOAT_EQ(sys.vz[1], -13.141091) 	<< "vz[1] is not correct";

}

