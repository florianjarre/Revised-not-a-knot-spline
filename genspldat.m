% genspldat.m
% GENerate some SPLine DATa x and f
% x and f must be column vectors
% f must be the same as in compspline if the latter is used !!!

n = 96; 
%a = pi/4;
%a = 0;
%b = a+pi;
a = -1;
b = 4;
a
b
n
%x = sort(rand(n,1)); disp('random knots')
x = 1:n; x = x.'; disp('equidistant knots')
xm = min(x); xM = max(x);
x = (x-xm)/(xM-xm);
x = a+(b-a)*x;     

% form some function 
%f = sin(x);   disp('sin')
%f = x.^4;     disp('x^4')
%f = x.^5;     disp('x^5')
%f = 1./(1+x.^2); disp('1./(1+x.^2)')
%f = exp(-1./(x.^2+1e-30)); disp('exp(-1./(x.^2))')
%f = atan(x); disp('atan')
f = 1./(exp(-x)+1); disp('logistic function')
%f = exp(sqrt(x)); disp('exp(sqrt(x))')