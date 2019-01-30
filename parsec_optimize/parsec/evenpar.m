%% parameter range
%{
rle -0.5~0.5
Xup 0.2~0.6
Yup 0~1
YXXup 0~-20
Xlow 0.2~0.6
Ylow 0~1
YXXlow 0~-20
yte -0.5~0.5
delta yte 0~0.1
alpha (tail pointing)
beta (tail crossing angle) 0~90

%}
%% 
scale = [1 1]
shift = [-0.5 0]
par = [0.1 0.45 0.155 -0.2 0.45 0.155 -0.2 0 0 0 10];
% par = [0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,-0.1039,0,-2.7791,9.2496];
points = parsecpoints(par)