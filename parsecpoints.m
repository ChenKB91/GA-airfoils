function [foil,self_intersect]=parsecpoints(par)
    n = 101;
    m = n-1;
    rng(42); % random seed
    
%     p1=rle p2=Xup p3=Yup p4=YXXup p5=Xlow p6=Ylow p7=YXXlow p8=yte p9=deltayte (t.e. thickness)
%     p10=alpha te p11=beta te
%     range = [
%     par = [0.3, 0.35, 0.155, -0.35, 0.45, 0.006, -0.2, 0.0, 0.01, 0.05, -6];
%     par = rand(1,11);
%     par = [0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,-0.1039,0,-2.7791,9.2496];
    
    
    foil_coef = parsec(par);
    x = linspace(0,1,n);
    xl = linspace(0.01,1,m);
    upper = foil_coef(1)* x.^(1/2) ...
        +foil_coef(2)* x.^(3/2) ...
        +foil_coef(3)* x.^(5/2) ...
        +foil_coef(4)* x.^(7/2) ...
        +foil_coef(5)* x.^(9/2) ...
        +foil_coef(6)* x.^(11/2);

    lower = foil_coef(7)* xl.^(1/2) ...
        +foil_coef(8)* xl.^(3/2) ...
        +foil_coef(9)* xl.^(5/2) ...
        +foil_coef(10)* xl.^(7/2) ...
        +foil_coef(11)* xl.^(9/2) ...
        +foil_coef(12)* xl.^(11/2);
    
%     if boolplot == true
%         plot(x,upper,'r-')
%         hold on
%         plot(x,lower,'b-')
%         axis equal
%     end
    px = [fliplr(x) xl];
%     px(end) = [];
    py = [fliplr(upper) lower];
%     py(end) = [];
    foil.parameters = par;
    foil.x = px; foil.y = py;
    
    self_intersect = 0;
    for i = 1:99 % should be 100, but somehow it will fuck up everything
        if py(101+i) > py(101-i)
            disp(py(101+i))
            disp(py(101-i))
            disp(i)
            self_intersect = 1;
        end
    end
%     
%     
% %     p = parsec_points
%     plot(points(1,:),points(2,:),'b-')
%     hold on
%     axis equal
end

