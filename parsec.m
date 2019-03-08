%% PARSEC
% Author: Sean Wu, UC Davis MAE Dept. 
% Date:   Mar. 2015

% CC BY-SA, Sean Wu 2015

%{
This function determines a=[a1, a2, ...an] to solve the airfoil polynomial.
Zn=an(p)*X^(n-1/2), where n is the number of coordinates for the upper or
lower surface. 

Input is a vector of PARSEC parameters p=[p1, p2, ...pn] where
p1=rle         
p2=Xup
p3=Yup
p4=YXXup
p5=Xlow
p6=Ylow
p7=YXXlow
p8=yte
p9=delta yte (t.e. thickness)
p10=alpha te
p11=beta te
%}

function a=parsec(p)

% Define matricies
c1=[1,1,1,1,1,1];
c2=[p(2)^(1/2),p(2)^(3/2),p(2)^(5/2),p(2)^(7/2),p(2)^(9/2),p(2)^(11/2)];
c3=[1/2, 3/2, 5/2, 7/2, 9/2, 11/2];
c4=[(1/2)*p(2)^(-1/2), (3/2)*p(2)^(1/2),(5/2)*p(2)^(3/2),(7/2)...
    *p(2)^(5/2),(9/2)*p(2)^(7/2),(11/2)*p(2)^(9/2)];
c5=[(-1/4)*p(2)^(-3/2),(3/4)*p(2)^(-1/2),(15/4)*p(2)^(1/2),(35/4)...
    *p(2)^(3/2),(63/4)*p(2)^(5/2),(99/4)*p(2)^(7/2)];
c6=[1,0,0,0,0,0];

Cup=[c1; c2; c3; c4; c5; c6];

c7=[1,1,1,1,1,1];
c8=[p(5)^(1/2),p(5)^(3/2),p(5)^(5/2),p(5)^(7/2),p(5)^(9/2),p(5)^(11/2)];
c9=[1/2, 3/2, 5/2, 7/2, 9/2, 11/2];
c10=[(1/2)*p(5)^(-1/2), (3/2)*p(5)^(1/2),(5/2)*p(5)^(3/2),(7/2)...
    *p(5)^(5/2),(9/2)*p(5)^(7/2),(11/2)*p(5)^(9/2)];
c11=[(-1/4)*p(5)^(-3/2),(3/4)*p(5)^(-1/2),(15/4)*p(5)^(1/2),(35/4)...
    *p(5)^(3/2),(63/4)*p(5)^(5/2),(99/4)*p(5)^(7/2)];
c12=[1,0,0,0,0,0];

Clo=[c7; c8; c9; c10; c11; c12];

bup=[p(8)+p(9)/2;p(3);tand(-p(10)+p(11)/2);0;p(4);(sqrt(2*p(1)))];
blo=[p(8)-p(9)/2;p(6);tand(p(10)+p(11)/2);0;p(7);-(sqrt(2*p(1)))];

% Solve system of equations: C x a=b
aup=linsolve(Cup,bup); 
alower=linsolve(Clo,blo);

% Format a
% a = zeros(,);
a(:,1)=aup;
a(7:12,1)=alower;

end

% Sean Wu, 2015