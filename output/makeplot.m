clc;clear all

load('RedBlue.mat')
% -------- data sctructure ---------------------------------------------
% data:
% data.parameters : it         : snapshot time steps
%                   m          : number of grid point in grid x-direction
%                   n          : number of grid point in grid y-direction
%                   mg         : number of grid levels
%                   nb         : number of immersed boundary points
%                   ncr        : number of centers of rotation
%                   output_sf  : output streamfunction or not
%                   rey        : Reynolds number
%                   dt         : time increament
%                   len        : length of the first domain in grid
%                                x-direction
%                   offsetx    : offset of grid in grid x-direction
%                   offsety    : offset of grid in grid y-direction
%                   delta      : grid spacing on 1st domain (= len/m)
%data.levels(mg) :  coord
%                   qn   : un, vn
%                   q0pn : un, vn
%
% -------- General Parameters ---------------------------------------------

% Number of the time step of the file: this identifies which file should be read
istart = 1000;
iend = 1000;
inv = 100;
%Change colormap
load('RedBlue.mat');
dir = pwd;
it = iend;
pressure = 1;


for iph = istart:inv:iend

    fprintf('ploting figure at iph = %d \n',iph)

    [data] = getdata(dir,iph,pressure);
    % Range for plots
    range = [-1 2 -1 1];

    % Contour maximum values and number of contour levels
    % Vorticity
    cmax_w = 5; %sphere: 10;
    clev_w = 51;
    cmax_s = 4;
    clev_s = 40;

    dt = data.parameters.dt;
    % vy = 0.3*sin(2*pi*0.2*dt*it)
    % yshift = 0.3/(2*pi*0.2)*(-cos(2*pi*0.2*iph*dt));
    yshift = 0;

    figure(1)
    subplot(2,1,1)
    % vorticity field with streamlines
    mg = data.parameters.mg;
    for k = mg:-1:1
        xn = data.levels(k).coord_lab.xn;
        yn = data.levels(k).coord_lab.yn+yshift;
        wn = data.levels(k).wn;
        wn_plot = wn;
        for i = 1:size(wn,1)
            for j = 1:size(wn,2)
                if wn(i,j) > cmax_w
                    wn_plot(i,j) = cmax_w;
                elseif  wn(i,j) < -cmax_w
                    wn_plot(i,j) = -cmax_w;
                end
            end
        end
        warning('off','all');
        contourf(xn,yn,wn_plot,-cmax_w:2*cmax_w/clev_w:cmax_w,'edgecolor','none');
        axis equal;colorbar;hold on;axis(range);set(gca,'CLim',[-cmax_w cmax_w]);
        colormap(cmap);
        % plotting streamlines
        if k==3
            sn = data.levels(k).sn;
            xn_temp = xn;
            yn_temp = yn;
        end
        if k==1
            contour(xn_temp,yn_temp,sn,-cmax_s:2*cmax_s/clev_s:cmax_s,'edgecolor','k');hold on;
        end
        time = iph*dt;
        title(['Vorticity Field at t = ' num2str(time,'%2.1f'),' '],'fontsize',20)
    end
    % fill the bodies
    xb = data.ibdata.xb_lab;
    yb = data.ibdata.yb_lab+yshift;
    codeb = data.ibdata.codeb;
    for i = min(codeb):max(codeb)
        fill(xb(codeb==i),yb(codeb==i),'k');hold on
    end

    % plot a circle
    %theta = 0:2*pi/100:2*pi; rad = 1.5;
    %plot(rad*cos(theta),rad*sin(theta),'--k');hold on

    set(gca,'fontsize',18)
    hold off
    %

    % Pressure
    % Range for plots
    %range = [-2 2 -2 2];

    cmax_p = 1;
    clev_p = 51;

    subplot(2,1,2)
    % pressure field with streamlines
    mg = data.parameters.mg;
    for k = mg:-1:1
        xn = data.levels(k).coord_lab.xn;
        yn = data.levels(k).coord_lab.yn+yshift;
        pn = data.levels(k).pn;
        pn_plot = pn;
        for i = 1:size(pn,1)
            for j = 1:size(pn,2)
                if pn(i,j) > cmax_p
                    pn_plot(i,j) = cmax_p;
                elseif  pn(i,j) < -cmax_p
                    pn_plot(i,j) = -cmax_p;
                end
            end
        end
        warning('off','all');
        contourf(xn,yn,pn_plot,-cmax_p:2*cmax_p/clev_p:cmax_p,'edgecolor','none');
        axis equal;colorbar;hold on;axis(range);set(gca,'CLim',[-cmax_p cmax_p]);
        colormap(cmap);
        %plotting streamlines
        if k==1
            contour(xn_temp,yn_temp,sn,-cmax_s:2*cmax_s/clev_s:cmax_s,'edgecolor','k');hold on;
        end
        title(['Pressure Field at t = ' num2str(time,'%2.1f'),' '],'fontsize',20)
    end
    % fill the bodies
    xb = data.ibdata.xb_lab;
    yb = data.ibdata.yb_lab+yshift;
    codeb = data.ibdata.codeb;
    for i = min(codeb):max(codeb)
        fill(xb(codeb==i),yb(codeb==i),'k');hold on
    end
    % plot a circle
    %theta = 0:2*pi/100:2*pi; rad = 1.5;
    %plot(rad*cos(theta),rad*sin(theta),'--k');hold on

    set(gca,'fontsize',18)
    hold off


%     subplot(2,2,2)
%     plot(ang_rate(:,1)*dt,mod(ang_rate(:,2)/(2*pi),1),'b');hold on
%     plot(iph*dt,mod(data.constraint(1)/(2*pi),1),'ob',...
%          'markersize',8,'markerfacecolor','b');hold on
%     axis([0 ang_rate(end,1)*dt 0 1])
%     xlabel('time','fontsize',18);
%     ylabel('angle','fontsize',18);
%     set(gca,'fontsize',18)
%     hold off
%
%     subplot(2,2,4)
%     plot(ang_rate(:,1)*dt,ang_rate(:,3),'b');hold on
%     plot(iph*dt,data.constraint(2),'ob',...
%          'markersize',8,'markerfacecolor','b');hold on
%     axis([0 ang_rate(end,1)*dt 0 0.6])
%     xlabel('time','fontsize',18);
%     ylabel('angular velocity','fontsize',18);
%     set(gca,'fontsize',18)
%     hold off

    %
    set(gcf,'color','w')
    saveas(gcf,['fig',num2str(iph/inv)],'png');

end

%%

fprintf('ploting lift & drag curve ...\n')

force = load([pwd '/responds/force.dat']);
dt = data.parameters.dt;
tt = force(:,1)*dt;
drag = force(:,4);
lift = force(:,5);
frange = [0 10 0 5];

figure(2)
subplot(1,2,1)
plot(tt,drag,'-r');hold off
grid on;axis(frange)
xlabel('Time','fontsize',20)
ylabel('Drag','fontsize',20)
set(gca,'fontsize',18)

subplot(1,2,2)
plot(tt,lift,'-b');hold off
grid on;axis(frange)
xlabel('Time','fontsize',20)
ylabel('Lift','fontsize',20)
set(gca,'fontsize',18)

set(gcf,'color','w')
disp('saving....')
saveas(gcf,'force','png');

fprintf('Done. \n')