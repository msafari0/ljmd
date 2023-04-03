//I added this
#include <gtest/gtest.h>
#include "../include/ljmd.h"
#include "../include/mdlib.h"
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//test forces.c for the function force, which is used to compute the forces on each atom of a system of 2 atoms
TEST(forces, force) {
    mdsys_t sys;
    sys.natoms = 2;
    sys.rx[0] = 0.0;
    sys.rx[1] = 0.0;
    sys.ry[0] = 0.0;
    sys.ry[1] = 0.0;
    sys.rz[0] = 0.0;
    sys.rz[1] = 0.0;
    sys.fx[0] = 0.0;
    sys.fx[1] = 0.0;
    sys.fy[0] = 0.0;
    sys.fy[1] = 0.0;
    sys.fz[0] = 0.0;
    sys.fz[1] = 0.0;
    sys.mass = 39.948;
    sys.epsilon = 0.2379;
    sys.sigma = 3.405;
    sys.rcut = 8.5;
    sys.box = 17.1580;
    force(&sys);
    EXPECT_EQ(0.0, sys.fx[0]);
    EXPECT_EQ(0.0, sys.fx[1]);
    EXPECT_EQ(0.0, sys.fy[0]);
    EXPECT_EQ(0.0, sys.fy[1]);
    EXPECT_EQ(0.0, sys.fz[0]);
    EXPECT_EQ(0.0, sys.fz[1]);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}