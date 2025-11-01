
% comparethree_rev.m
% compare the natural spline, the nak spline, the Q spline 
% and the RNAK spline

clear
clc
genspldat
hold off

if 0 > -1
natspline
compspline % must be the same function as in genspldat
plot(t,ft-s,':') % dotted
dist_to_f_nat = dist_to_f
hold on
end

nakspline
compspline % must be the same function as in genspldat
plot(t,ft-s,'-') % solid
dist_to_f_nak = dist_to_f
hold on

if 0 > -1
qspline
compspline % must be the same function as in genspldat
plot(t,ft-s,'-.') % dashdot
dist_to_f_qspline = dist_to_f
hold on
end

rnakspline
compspline % must be the same function as in genspldat
plot(t,ft-s,'--') % dashed
dist_to_f_RNAK = dist_to_f
hold on

   
