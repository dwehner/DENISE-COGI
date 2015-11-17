#!/bin/csh

set x=1
set r

while ( $x < 101)

set r=`nawk 'BEGIN {srand();print int(rand()*100000)}'`

suaddnoise sn=200 seed=$r < su/MARMOUSI_spike/DENISE_MARMOUSI_y.su.shot$x f=3,15,15,30 amps=0,1,1,0 > su/MARMOUSI_spike_noise/DENISE_MARMOUSI_y.su.shot$x
suaddnoise sn=200 seed=$r < su/MARMOUSI_spike/DENISE_MARMOUSI_x.su.shot$x f=3,15,15,30 amps=0,1,1,0 > su/MARMOUSI_spike_noise/DENISE_MARMOUSI_x.su.shot$x

set x = `expr $x + 1`

end


