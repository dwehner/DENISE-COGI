#!/bin/csh

set x=1
set folder='MARMOUSI_spike_BW_5Hz'

while ( $x < 101)

sufft < su/$folder/DENISE_MARMOUSI_y.su.shot$x | suamp mode=amp > su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_y.su.shot$x
sufft < su/$folder/DENISE_MARMOUSI_x.su.shot$x | suamp mode=amp > su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_x.su.shot$x

set x = `expr $x + 1`

end


susum su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_y.su.shot1 su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_y.su.shot2 > su/AMP-Spec_$folder/AMP_Spec_y_sum_test1

susum su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_x.su.shot1 su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_x.su.shot2 > su/AMP-Spec_$folder/AMP_Spec_x_sum_test1


set n=3
set m=1
set k=2

while ( $n < 101)

susum su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_y.su.shot$n su/AMP-Spec_$folder/AMP_Spec_y_sum_test$m > su/AMP-Spec_$folder/AMP_Spec_y_sum_test$k

susum su/AMP-Spec_$folder/AMP_Spec_DENISE_MARMOUSI_x.su.shot$n su/AMP-Spec_$folder/AMP_Spec_x_sum_test$m > su/AMP-Spec_$folder/AMP_Spec_x_sum_test$k


set n = `expr $n + 1`
set m = `expr $m + 1`
set k = `expr $k + 1`

end

mv su/AMP-Spec_$folder/AMP_Spec_y_sum_test99 su/AMP-Spec_$folder/AMP_Spec_y_all_shots
mv su/AMP-Spec_$folder/AMP_Spec_x_sum_test99 su/AMP-Spec_$folder/AMP_Spec_x_all_shots

rm su/AMP-Spec_$folder/AMP_Spec_y_sum_test*
rm su/AMP-Spec_$folder/AMP_Spec_x_sum_test*

sustack < su/AMP-Spec_$folder/AMP_Spec_y_all_shots > su/AMP-Spec_$folder/AMP_Spec_y_sum

sustack < su/AMP-Spec_$folder/AMP_Spec_x_all_shots > su/AMP-Spec_$folder/AMP_Spec_x_sum
