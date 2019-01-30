
function [data] = getdata(dir,it,pressure)

% Modified by Thibault Flinois on 23 Jan 2013

% Reads the data file at istep=it, and returns:
%    xb,yb :     The x,y coordinates of each body/actuator point
%    codeb :     A number identifying the body to which belongs each body/actuator point
%    xn,yn :     Coordinates of each grid point (Cartesian grid)
%    un,vn :     Velocity components at each grid point (based on fluxes from code, which are evaluated at cell faces)
%    wn :        Vorticity at grid point (based on circulation from code, which is evaluated at cell vertices)
%    sn :        Streamfunction at each grid point
%    pn :        Pressure at each grid point
%    filename :  Full name and path of the file read
%
% Input parameters:
%    dir :       Directory where the file is located
%    it :        Number of the time step of the file: this identifies which file should be read
%    lev :       Grid level where data is required
%    pressure :  1 if pressure data is required, 0 if not
%
% Note: The "*.var" files read here have been written to Fortran's
% unformatted format so in order for these (binary) files to be read, one
% needs to read several Fortran markers, which are of no interest here.
%
% -------------------------------------------------------------------------

% ---- Open file ----------------------------------------------------------
filename = [dir '/snapshots/ib' num2str(it,'%7.7i') '.var'];
%filename
fid = fopen(filename,'r');

% ---- Read first line of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    it = fread(fid,1,'int');       % time steps
    m = fread(fid,1,'int');           % cells in x direction
    n = fread(fid,1,'int');           % cells in y direction
    mg = fread(fid,1,'int');          % total number of grid levels
    nb = fread(fid,1,'int');          % total number of body points
    ncr = fread(fid,1,'int');         % total number of centers of rotation
    output_sf = fread(fid,1,'int' );  % output stream-function or not : 1: yes; 0: no
temp = fread(fid, 1, 'float32');      % Fortran marker

    data.parameters.it = it;       % time steps
    data.parameters.m = m;           % cells in x direction
    data.parameters.n = n;           % cells in y direction
    data.parameters.mg = mg;          % total number of grid levels
    data.parameters.nb = nb;          % total number of body points
    data.parameters.ncr = ncr;         % total number of centers of rotation
    data.parameters.output_sf = output_sf;  % output stream-function or not : 1: yes; 0: no


% ---- Read second line of file -------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    rey = fread(fid,1,'real*8' );     % Reynolds number
    dt = fread(fid,1,'real*8' );      % Size of time step
    len = fread(fid,1,'real*8' );     % Size of the smallest grid in the x direction
% Note: 'len' sets the size in the y direction too since the grid spacing is uniform
    offsetx = fread(fid,1,'real*8' ); % x-distance from lower-right corner smallest grid to origin
    offsety = fread(fid,1,'real*8' ); % y-distance of lower-right corner of smallest grid to origin
temp = fread(fid, 1, 'float32');      % Fortran marker
    data.parameters.rey = rey;     % Reynolds number
    data.parameters.dt = dt;      % Size of time step
    data.parameters.len = len;     % Size of the smallest grid in the x direction
% Note: 'len' sets the size in the y direction too since the grid spacing is uniform
    data.parameters.offsetx = offsetx; % x-distance from lower-right corner smallest grid to origin
    data.parameters.offsety = offsety; % y-distance of lower-right corner of smallest grid to origin
    data.parameters.delta = len./m;
% ---- Compute grid related parameters ------------------------------------

% If we are not considering the smallest grid level, we need to compute the
% grid spacing and x and y offsets for the current grid level.
%
% Note:  the cells in each grid level are twice as large in both directions
% as the ones of the previous grid level
for k = 1:mg
    fac = 2^(k-1);
    delta = len ./ m *fac;                                  % Grid spacing in both directions for current grid
    offx = 2^(k-1) * len/2 - len/2 + offsetx;             % Offset in x direction for current grid
    offy = 2^(k-1) * (n*len/m)/2 - (n*len/m)/2 + offsety; % Offset in y direction for current grid

    temp = fread(fid, 1, 'float32');      % Fortran marker
    q = fread(fid, (m+1)*(n)+(n+1)*(m), 'real*8');                              % Read fluxes at cell faces
    q0p = fread(fid, (m+1)*(n)+(n+1)*(m), 'real*8');                              % Read fluxes at cell faces
    nl_old = fread(fid, (m+1)*(n)+(n+1)*(m), 'real*8');                              % Read fluxes at cell faces
    temp = fread(fid, 1, 'float32');      % Fortran marker

    qx(:,:) = transpose(reshape(q(1:(m+1)*(n)),m+1,n))/delta;                       % Put x-comp of fluxes in array
    qy(:,:) = transpose(reshape(q((m+1)*(n)+1:(m+1)*(n)+(n+1)*(m)),m,n+1))/delta;   % Put y-comp of fluxes in array
    q0px(:,:) = transpose(reshape(q0p(1:(m+1)*(n)),m+1,n))/delta;                       % Put x-comp of fluxes in array
    q0py(:,:) = transpose(reshape(q0p((m+1)*(n)+1:(m+1)*(n)+(n+1)*(m)),m,n+1))/delta;   % Put y-comp of fluxes in array
    nl_old_x(:,:) = transpose(reshape(nl_old(1:(m+1)*(n)),m+1,n))/(delta^3);                       % Put x-comp of fluxes in array
    nl_old_y(:,:) = transpose(reshape(nl_old((m+1)*(n)+1:(m+1)*(n)+(n+1)*(m)),m,n+1))/(delta^3);   % Put y-comp of fluxes in array

    q0r = zeros((m+1)*(n)+(n+1)*(m),ncr);
    q0rx = zeros(n,m+1,ncr);
    q0ry = zeros(n+1,m,ncr);
    for i = 1:ncr
        temp = fread(fid, 1, 'float32');      % Fortran marker
        q0r(:,i) = fread(fid, (m+1)*(n)+(n+1)*(m), 'real*8');                              % Read fluxes at cell faces
        temp = fread(fid, 1, 'float32');      % Fortran marker
        q0rx(:,:,i) = transpose(reshape(q0r(1:(m+1)*(n),i),m+1,n))/delta;                       % Put x-comp of fluxes in array
        q0ry(:,:,i) = transpose(reshape(q0r((m+1)*(n)+1:(m+1)*(n)+(n+1)*(m),i),m,n+1))/delta;   % Put y-comp of fluxes in array
    end

    temp = fread(fid, 1, 'float32');      % Fortran marker
    omega(:,:) = transpose(reshape(fread(fid, (m-1)*(n-1), 'real*8'), m-1,n-1));
    temp = fread(fid, 1, 'float32');      % Fortran marker
    omega = omega / delta.^2;

    if (output_sf == 1)
        temp = fread(fid, 1, 'float32');      % Fortran marker
        stfn(:,:) = transpose(reshape(fread(fid, (m-1)*(n-1), 'real*8'), m-1,n-1));
        temp = fread(fid, 1, 'float32');      % Fortran marker

        % defined psi0', which psi0'(1,1)= 0
        stfn0 = zeros(n-1,m-1);
        for i = 2:m-1
            stfn0(1,i) = stfn0(1,i-1) - delta*q0py(2,i);
        end
        for j = 2:n-1
            stfn0(j,1:m-1) = stfn0(j-1,1:m-1) + delta*q0px(j,2:m);
        end

    end

    % ---- Interpolate all variables to the same grid -------------------------

    % Note: All variables are interpolated to the cell vertices

    % Create the grid that will be used for all variables
    [xn, yn] = meshgrid(delta:delta:(m-1)*delta, delta:delta:(n-1)*delta);
    xn = xn - offx;
    yn = yn - offy;

    % Grid for x-velocities (vertical cell faces)
    [xu, yu] = meshgrid(0:delta:m*delta, delta/2:delta:(n-0.5)*delta);
    xu = xu - offx;
    yu = yu - offy;

    % Grid for y-velocities (horizontal cell faces)
    [xv, yv] = meshgrid(delta/2:delta:(m-0.5)*delta, 0:delta:n*delta);
    xv = xv - offx;
    yv = yv - offy;

    % Grid for vorticity and streamfunction (cell vertices)
    [xw, yw] = meshgrid(delta:delta:(m-1)*delta, delta:delta:(n-1)*delta);
    xw = xw - offx;
    yw = yw - offy;

    % Interpolate all variables accordingly to xn, yn
    un = interp2(xu,yu,qx,xn,yn);
    vn = interp2(xv,yv,qy,xn,yn);
    u0pn = interp2(xu,yu,q0px,xn,yn);
    v0pn = interp2(xv,yv,q0py,xn,yn);
    for i=1:ncr
        u0rn(:,:,i) = interp2(xu,yu,q0rx(:,:,i),xn,yn);
        v0rn(:,:,i) = interp2(xv,yv,q0ry(:,:,i),xn,yn);
    end
    nl_un = interp2(xu,yu,nl_old_x,xn,yn);
    nl_vn = interp2(xv,yv,nl_old_y,xn,yn);

    data.levels(k).coord.xn = xn;
    data.levels(k).coord.yn = yn;
    data.levels(k).qn.un = un;
    data.levels(k).qn.vn = vn;
    data.levels(k).q0pn.un = u0pn;
    data.levels(k).q0pn.vn = v0pn;
    data.levels(k).nonlinear.un = nl_un;
    data.levels(k).nonlinear.vn = nl_vn;

    u0rn = zeros(m-1,n-1,ncr);
    v0rn = zeros(m-1,n-1,ncr);
    for i = 1:ncr
        u0rn(:,:,i) = interp2(xu,yu,q0rx(:,:,i),xn,yn);
        v0pn(:,:,i) = interp2(xv,yv,q0ry(:,:,i),xn,yn);
        data.levels(k).q0rn(i).un = u0rn(:,:,i);
        data.levels(k).q0rn(i).vn = v0rn(:,:,i);
    end
    data.levels(k).wn = interp2(xw,yw,omega,xn,yn);
    if (output_sf == 1)
        disp('hi')
        [xw0, yw0] = meshgrid(-delta:delta:delta,-delta:delta:delta);
        s0 = interp2(xw,yw,stfn0,xw0,yw0);
        sn = interp2(xw,yw,stfn+stfn0,xn,yn);
        data.levels(k).sn = sn-s0(2,2);
    end

end

% ---- Read third line of file --------------------------------------------
if ( nb > 0)
% --> Body coordinates and body velocity relative to the grid
    temp = fread(fid, 1, 'float32');      % Fortran marker
    xyb = fread(fid, 2*nb, 'real*8');   % Read x and y coordinates of each body point
    xb = xyb(1:nb);                     % Create x coordinates vector for body points
    yb = xyb(nb+1:2*nb);                % Create y coordinates vector for body points
    data.ibdata.xb = xb;
    data.ibdata.yb = yb;

    vxyb = fread(fid, 2*nb, 'real*8');    % Read x and y components of velocity of each body point
    vxb = vxyb(1:nb);                     % Create x-velocity for body points
    vyb = vxyb(nb+1:2*nb);                % Create y-velocity for body points
    data.ibdata.vxb = vxb;
    data.ibdata.vyb = vyb;

    fxyb = 2*fread(fid, 2*nb, 'real*8')/dt;% Read x and y components of force coeff. of each body point
    fxb = fxyb(1:nb);                     % Create x-force coeff. for body points
    fyb = fxyb(nb+1:2*nb);                % Create y-force coeff.for body points
    data.ibdata.fxb = fxb;
    data.ibdata.fyb = fyb;
    temp = fread(fid, 1, 'float32');

    temp = fread(fid, 1, 'float32');
    codeb = fread(fid, nb, 'int');   % Read body number corresponding to each body point
    data.ibdata.codeb = codeb;
    temp = fread(fid, 1, 'float32');      % Fortran marker
end
% ---- Read fifth line of file -------------------------------------------

    temp = fread(fid, 1, 'float32');
    xcm = fread(fid,2,'real*8' );
    vcm = fread(fid,2,'real*8' );
    theta_g2l = fread(fid,1,'real*8' );
    omega_g2l = fread(fid,1,'real*8' );
    temp = fread(fid, 1, 'float32');      % Fortran marker

    for i = 1:ncr
        temp = fread(fid, 1, 'float32');      % Fortran marker
        data.cm_motion.rot_cm(i).rox = fread(fid,1,'real*8' );   % rotating angle of the frame
        data.cm_motion.rot_cm(i).roy = fread(fid,1,'real*8' );   % rotating angle of the frame
        data.cm_motion.rot_cm(i).theta = fread(fid,1,'real*8' );   % rotating angle of the frame
        data.cm_motion.rot_cm(i).omega = fread(fid,1,'real*8' );   % rotating angle of the frame
        temp = fread(fid, 1, 'float32');      % Fortran marker
    end

    data.cm_motion.trans_cm.xcm = xcm;
    data.cm_motion.trans_cm.vcm = vcm;
    data.cm_motion.theta_g2l = theta_g2l;   % rotating angle of the frame
    data.cm_motion.omega_g2l = omega_g2l;


    % ---- Rotate the grid and the bodies to the lab-frame --------------------
    xcm_lab = xcm;
    xcm_lab(1) = xcm(1)*cos(theta_g2l) - xcm(2)*sin(theta_g2l) ;
    xcm_lab(2) = xcm(1)*sin(theta_g2l) + xcm(2)*cos(theta_g2l) ;
    data.cm_motion.trans_cm.xcm_lab = xcm_lab;

    vcm_lab = vcm;
    vcm_lab(1) = vcm(1)*cos(theta_g2l) - vcm(2)*sin(theta_g2l) ;
    vcm_lab(2) = vcm(1)*sin(theta_g2l) + vcm(2)*cos(theta_g2l) ;
    data.cm_motion.trans_cm.vcm_lab = vcm_lab;

    for k = 1:mg
        un = data.levels(k).qn.un;
        vn = data.levels(k).qn.vn;
        un_lab = un*cos(theta_g2l) - vn*sin(theta_g2l) ;
        vn_lab = un*sin(theta_g2l) + vn*cos(theta_g2l) ;
        data.levels(k).qn_lab.un = un_lab;
        data.levels(k).qn_lab.vn = vn_lab;

        u0pn = data.levels(k).q0pn.un;
        v0pn = data.levels(k).q0pn.vn;
        u0pn_lab = u0pn*cos(theta_g2l) - v0pn*sin(theta_g2l) ;
        v0pn_lab = u0pn*sin(theta_g2l) + v0pn*cos(theta_g2l) ;
        data.levels(k).q0pn_lab.un = u0pn_lab;
        data.levels(k).q0pn_lab.vn = v0pn_lab;

        xn = data.levels(k).coord.xn;
        yn = data.levels(k).coord.yn;

        % user-define rotation
        for i = 1:ncr

            rox = data.cm_motion.rot_cm(i).rox;
            roy = data.cm_motion.rot_cm(i).roy;
            %theta = data.cm_motion.rot_cm(i).theta;

            u0rn = data.levels(k).q0rn(i).un;
            v0rn = data.levels(k).q0rn(i).vn;
            u0rn_lab = u0rn*cos(theta_g2l) - v0rn*sin(theta_g2l) ;
            v0rn_lab = u0rn*sin(theta_g2l) + v0rn*cos(theta_g2l) ;
            data.levels(k).q0rn_lab(i).un = u0rn_lab;
            data.levels(k).q0rn_lab(i).vn = v0rn_lab;

            xnp = rox + (xn-rox)*cos(theta_g2l) - (yn-roy)*sin(theta_g2l);
            ynp = roy + (xn-rox)*sin(theta_g2l) + (yn-roy)*cos(theta_g2l);

            data.levels(k).coord_lab.xn = xnp;
            data.levels(k).coord_lab.yn = ynp;

        end

    end

    for i = 1:ncr

        rox = data.cm_motion.rot_cm(i).rox;
        roy = data.cm_motion.rot_cm(i).roy;

        xbp = rox + (xb-rox)*cos(theta_g2l) - (yb-roy)*sin(theta_g2l);
        ybp = roy + (xb-rox)*sin(theta_g2l) + (yb-roy)*cos(theta_g2l);

        data.ibdata.xb_lab = xbp;
        data.ibdata.yb_lab = ybp;

    end

% ---- Done reading file --------------------------------------------------
fclose(fid);


% ---- Pressure Computations ----------------------------------------------
if(pressure == 1)
    % Open file
    filename = [dir '/snapshots/pressure' num2str(it,'%7.7i') '.var'];
    fid = fopen(filename,'r');
% --> Pressure
    % Read pressure at grid cell centres, for each grid level, and output in m-by-n arrays
    for k = 1:mg
        temp = fread(fid, 1, 'float32');      % Fortran marker
        p(:,:) = transpose(reshape(fread(fid, (m)*(n), 'real*8'), m,n));
        temp = fread(fid, 1, 'float32');      % Fortran marker

        temp = fread(fid, 1, 'float32');      % Fortran marker
        pdyn(:,:) = transpose(reshape(fread(fid, (m)*(n), 'real*8'), m,n));
        temp = fread(fid, 1, 'float32');      % Fortran marker

        fac = 2^(k-1);
        delta = len ./ m *fac;                                  % Grid spacing in both directions for current grid
        offx = 2^(k-1) * len/2 - len/2 + offsetx;             % Offset in x direction for current grid
        offy = 2^(k-1) * (n*len/m)/2 - (n*len/m)/2 + offsety; % Offset in y direction for current grid

        % Create the grid that will be used for all variables
        [xn, yn] = meshgrid(delta:delta:(m-1)*delta, delta:delta:(n-1)*delta);
        xn = xn - offx;
        yn = yn - offy;

        % Grid for Pressure (cell centres)
        [xp yp]=meshgrid(delta:delta:m*delta, delta:delta:n*delta);
        xp = xp - offx;
        yp = yp - offy;

        pn = interp2(xp,yp,p,xn,yn);
        pn_dyn = interp2(xp,yp,pdyn,xn,yn);
        data.levels(k).pn = 2*pn;
        data.levels(k).pn_dyn = 2*pn_dyn;
    end
    % Done reading
	fclose(fid);

end

end
