function [newpt, self_intersect] = evenpar(par)
%% parameter range
%{
rle 0~0.5
Xup 0.2~0.6
Yup 0~1
YXXup 0~-20
Xlow 0.2~0.6
Ylow 0~1
YXXlow 0~-20
yte -0.5~0.5
delta yte   0   % for the fortran code be able to run
alpha (tail pointing) 45~-45
beta (tail crossing angle) 0~90

%}
%%
% min = [0   0.3 0   -20 0.2 -0.3 -20 -0.5 0 -45 0];
% max = [0.1 0.6 0.3 0   0.6  0     0 0.5 0 45 90];
min = [0.02  0.32 0.077 -0.65 0.15 0.02 -0.60 -0.1 0 -4.55 15];
max = [0.023 0.37 0.080 -0.63 0.19 0.05 -0.75  0.1 0 -4.90 15.1];
% par = [0.015 0.45 0.07 -1.75 0.3 -0.06 0.45 0 0 13.63 32.16];
% par = [0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,-0.1039,0,-2.7791,9.2496];
[points, self_intersect] = parsecpoints(par);

pt = interparc(100, points.x, points.y, 'csape');
% pt(100,:) = []
targetd = 0.015;
d = pdist([pt(50,:);pt(51,:)]);
newn = 100;
while 1
    newn = round(newn*d/targetd);
    %             disp(newn);
    newpt = interparc(newn, points.x, points.y);
    newpt(newn,:) = [];
    shift = 0;
    newpt(:,1) = newpt(:,1)-shift;
    
    %             plot(newpt(:,1),newpt(:,2),'bo')
    %             hold on
    %             axis equal
    d = pdist([newpt(1,:);newpt(2,:)]);
    if abs((d - targetd)/targetd) < 0.01
        break
    else
%         fprintf(logfile, '  d = %f\n, regen', d);
    end
end
%
[~,i] = maxk(newpt(:,1), 4);
p1 = newpt(i(1),:);
p2 = newpt(i(2),:);
p3 = newpt(i(3),:);
p4 = newpt(i(4),:);
state = [p2(2) p3(2) p4(2)]>p1(2);
excp = isequal(state,[1 0 0]) || isequal(state,[0 0 1]) || isequal(state,[1 0 1]) || isequal(state,[0 1 0]);

if excp
    pa = p2;
    pb = p3;
else
    pa = p3;
    pb = p4;
end


vec_vert_len = sqrt(0.015^2 - (norm(pa-pb)/2)^2);
vec_vert = vec_vert_len * fliplr((pa-pb)/norm(pa-pb)) .* [-1, 1];
if pa(2) > pb(2)
    vec_vert = -vec_vert;
end
tmppt = (pa+pb)/2 + vec_vert;

if excp
    newpt(i(1),:) = tmppt;
else
    newpt(i(1),:) = [];
    newpt(i(2),:) = tmppt;
end
%}
% plot(newpt(:,1),newpt(:,2),'bo')
% hold on
% axis equal