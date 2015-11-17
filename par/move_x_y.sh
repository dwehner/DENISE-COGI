#!/bin/csh

set x=1
while ( $x < 101)

mv su/DENISE_MARMOUSI_x.su.shot$x.it1 su/TEST_steep_fault/DENISE_FAULT_x.su.shot$x
mv su/DENISE_MARMOUSI_y.su.shot$x.it1 su/TEST_steep_fault/DENISE_FAULT_y.su.shot$x

set x = `expr $x + 1`

end
