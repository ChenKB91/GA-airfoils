function points=parsecpoints(par)
    n = 100;
    rng(42); % random seed
    
%     p1=rle p2=Xup p3=Yup p4=YXXup p5=Xlow p6=Ylow p7=YXXlow p8=yte p9=deltayte (t.e. thickness)
%     p10=alpha te p11=beta te
%     range = [
%     par = [0.3, 0.35, 0.155, -0.35, 0.45, 0.006, -0.2, 0.0, 0.01, 0.05, -6];
%     par = rand(1,11);
%     par = [0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,-0.1039,0,-2.7791,9.2496];
    
    
    foil_coef = parsec(par);
    x = linspace(0,1,n);
    upper = foil_coef(1)* x.^(1/2) ...
        +foil_coef(2)* x.^(3/2) ...
        +foil_coef(3)* x.^(5/2) ...
        +foil_coef(4)* x.^(7/2) ...
        +foil_coef(5)* x.^(9/2) ...
        +foil_coef(6)* x.^(11/2);

    lower = foil_coef(7)* x.^(1/2) ...
        +foil_coef(8)* x.^(3/2) ...
        +foil_coef(9)* x.^(5/2) ...
        +foil_coef(10)* x.^(7/2) ...
        +foil_coef(11)* x.^(9/2) ...
        +foil_coef(12)* x.^(11/2);

    plot(x,upper,'bo-')
%     hold on
%     plot(x,-lower,'ro-')
%     axis equal

    px = [x fliplr(x)];
    py = [upper fliplr(-lower)];
    points = [px;py];
    
    
%     p = parsec_points
    plot(points(1,:),points(2,:),'b-')
    hold on
    axis equal
end

