% qspline.m
% Compute a Q-spline for approximating a function f just based
% on its function values at given points. 
% Approximate f''(x(1)) via cubic interpolation and improve f'' using 
% the fourth order correction. Bound the correction by a safeguard. 
%
% Input:  column vectors x, f of same dimension n with n >= 6 (!!!)
%         x must be strictly increasing, x(i+1) > x(i) for all i.
% Output: (n-1) x 4 - matrix S defining the Q-spline.
%         for x in [x(i),x(i+1)], the spline value s(x) is   (Horner)
%         S(i,1)+ (x-x(i))*(S(i,2)+ (x-x(i))*(S(i,3)+ (x-x(i))*(S(i,4))))
%
% for a test, call genspldat.m first to generate some test data

% Florian Jarre, Oct 29, 2025
% For testing only, no guarantee for correctness

% set up the tri-diagonal system for the moments
n      = length(x);
h      = x(2:n)-x(1:n-1);   % h(j)  = x(j+1)-x(j)  (in R^{n-1})
hp     = h(1:n-2)+h(2:n-1); % hp(j) = x(j+2)-x(j)  (in R^{n-2})
mu     = h(1:n-2)./hp; %                           (in R^{n-2})
lambda = h(2:n-1)./hp; %                           (in R^{n-2})
A = 2*speye(n-2)+spdiags([mu(2:n-2);    0],-1,n-2,n-2)+...
                 spdiags([0;lambda(1:n-3)], 1,n-2,n-2);
% elements mu(n-1) or lambda(0) are not used....
% same as
% A      = 2*eye(n-2)+diag(mu(2:n-2),-1)+diag(lambda(1:n-3),1);

% form table of divided differences explicitly
% (overkill, only the end points f4(1) and f4(n-4) are neded)
f0 = f;
f1 = (f0(2:n)  -f0(1:n-1))./(x(2:n)-x(1:n-1)); % n-1  elements
f2 = (f1(2:n-1)-f1(1:n-2))./(x(3:n)-x(1:n-2)); % n-2  elements
f3 = (f2(2:n-2)-f2(1:n-3))./(x(4:n)-x(1:n-3)); % n-3  elements
f4 = (f3(2:n-3)-f3(1:n-4))./(x(5:n)-x(1:n-4)); % n-4  elements
 

% compute natural spline S1 and two zero-interpolating splines S2, S3
rhs1 = 6*f2;                                   % natural spline
rhs2 = zeros(n-2,1); rhs2(1)   = -mu(1);       % (moment M1 = 1)
rhs3 = zeros(n-2,1); rhs3(n-2) = -lambda(n-2); % (moment Mn = 1)

% compute moments of S1, S2, S3
M1 = A \ rhs1;  M1 = [0;M1;0];
M2 = A \ rhs2;  M2 = [1;M2;0];
M3 = A \ rhs3;  M3 = [0;M3;1];

% for j = 1:n-1 we have for the interval [x(j),x(j+1)] of length h(j)
% c(j) = f1(j) - (M(j+1)-M(j))*h(j)/6
% d(j) = f0(j) -  M(j)*h(j)^2/6;
% s(x) = (M(j+1)*(x-x(j))^3 + M(j)*(x(j+1)-x)^3)/(6*h(j))
%        +c(j)*(x-x(j)) + d(j)
% = ( (M(j+1)-M(j)) /(6*h(j)) ) *(x-x(j))^3
%   + 0.5*M(j)                  *(x-x(j))^2
%   + (c(j) - 0.5*M(j)*h(j))    *(x-x(j))
%   + d(j) + (M(j)*h(j)^2) /6

% compute representations of s_{nat}=S1, S2, and S3
c1 = f1 - (M1(2:n)-M1(1:n-1)).*h/6;
d1 = f0(1:n-1) - (M1(1:n-1).*h.^2)/6;
S1 = zeros(n-1,4);
S1(:,1) = d1 + (M1(1:n-1).*h.^2) /6; % this is f0 !!
S1(:,2) = c1 - 0.5*M1(1:n-1).*h;
S1(:,3) = 0.5*M1(1:n-1);
S1(:,4) = (M1(2:n)-M1(1:n-1)) ./(6*h);

c2 =  - (M2(2:n)-M2(1:n-1)).*h/6;
d2 =  - (M2(1:n-1).*h.^2)/6;
S2 = zeros(n-1,4);
S2(:,1) = d2 + (M2(1:n-1).*h.^2) /6; % this is zero !!
S2(:,2) = c2 - 0.5*M2(1:n-1).*h;
S2(:,3) = 0.5*M2(1:n-1);
S2(:,4) = (M2(2:n)-M2(1:n-1)) ./(6*h);

c3 =  - (M3(2:n)-M3(1:n-1)).*h/6;
d3 =  - (M3(1:n-1).*h.^2)/6;
S3 = zeros(n-1,4);
S3(:,1) = d3 + (M3(1:n-1).*h.^2) /6; % this is zero !!
S3(:,2) = c3 - 0.5*M3(1:n-1).*h;
S3(:,3) = 0.5*M3(1:n-1);
S3(:,4) = (M3(2:n)-M3(1:n-1)) ./(6*h);


f5l = (f4(2)-f4(1))/(x(6)-x(1));  
% plain divided difference, f5l \approx f^{(5)}((x(3)+x(4))/2) / 120
f5r = (f4(n-4)-f4(n-5))./(x(n)-x(n-5));
% same for the other end point


% find estimates f''(x(1))=fppl and f''(x(n))=fppr   
% SAFEGURAD: possibly reduce the fourth order correction terms using
% deltarhol and deltarhol before adding the correction terms to f
if f4(1)*f5l > 0
    deltarhol = 2.5*(x(3)-x(2))*f5l/f4(1); % this is positive
else
    deltarhol = 0; % no adjustment 
end
%deltarhol = 0; % pure Q-spline without adjustment
rhol = f4(1)*max(0,1-deltarhol); % damping based on fifth derivative
% f4(1)\approx f^{(4)}(x(3))/24,   
% rhol \approx f^{(4)}(x(2)+x(3))/2) / 24 

if f4(n-4)*f5r < 0
    deltarhor = 2.5*(x(n-1)-x(n-2))*f5r/f4(n-4); % this is negative
else
    deltarhor = 0;
end
%deltarhor = 0; % pure Q-spline without adjustment
rhor = f4(n-4)*max(0,1+deltarhor); % damping based on fifth derivative


% form fhat and its table of divided differences for left end point
% cubic function on [x(1),x(4)] 
% for further testing, also values wrt. x(5), x(6) are computed 
fh0 = f0(1:6)  - rhol*(x(1:6)-x(1)).^4;     % function values of fhat
fh1 = (fh0(2:6)-fh0(1:5))./(x(2:6)-x(1:5)); % five elements
fh2 = (fh1(2:5)-fh1(1:4))./(x(3:6)-x(1:4)); % four  elements
fh3 = (fh2(2:4)-fh2(1:3))./(x(4:6)-x(1:3)); % three elements

% fh(x) \approx  fh0(1) + fh1(1)*(x-x(1)) 
%                       + fh2(1)*(x-x(1))*(x-x(2))  
%                       + fh3(1)*(x-x(1))*(x-x(2))*(x-x(3))

% cubic estimate of f''(x(1)) via fhat 
fppl = 2*fh2(1) + 2*fh3(1)*(2*x(1)-x(2)-x(3));


% form fhat and its table of divided differences for right end point
% cubic function on [x(n-3),x(n)]
% for further testing, also values wrt. x(n-5), x(n-4) are computed 
fh0 = f0(n-5:n) - rhor*(x(n-5:n)-x(n)).^4;        % function values fhat, 
fh1 = (fh0(2:6)-fh0(1:5))./(x(n-4:n)-x(n-5:n-1)); % five  elements
fh2 = (fh1(2:5)-fh1(1:4))./(x(n-3:n)-x(n-5:n-2)); % four  elements
fh3 = (fh2(2:4)-fh2(1:3))./(x(n-2:n)-x(n-5:n-3)); % three elements

% f(x) \approx  f0(3) + f1(3)*(x-x(n-3)) 
%                     + f2(3)*(x-x(n-3))*(x-x(n-2))  
%                     + f3(3)*(x-x(n-3))*(x-x(n-2))*(x-x(n-1))
% fpp(x6) \approx 2*fh2(3) + fh3(3) * (6*x(n)-2*(x(n-3)+x(n-2)+x(n-1)))

% cubic estimate of f''(x(6)) via fhat 
fppr = 2*fh2(3) + fh3(3)*(6*x(n)-2*(x(n-3)+x(n-2)+x(n-1)));

%fppl = 0; fppr = 0; % this leads to the natural spline



% Set up the coefficient system for the RQ-spline
A = [...
     2*S2(1,3),                      2*S3(1,3);...
     2*S2(n-1,3)+6*S2(n-1,4)*h(n-1), 2*S3(n-1,3)+6*S3(n-1,4)*h(n-1)...
     ];

rhs = [...
      fppl-2*S1(1,3);...
      fppr-2*S1(n-1,3)-6*S1(n-1,4)*h(n-1)...
      ];
              
z = A \ rhs;

S = S1+z(1)*S2+z(2)*S3;


