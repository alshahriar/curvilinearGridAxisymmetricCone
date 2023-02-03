%% 3D Semi-body conformal curvilinear structured grid generator for axisymmetric body
% Grid points are structured as ndgrid() (not meshgrid)
% But in curvilinear coordinates (not Cartesian grid)

%% Outputs
% Full grid in tecplot format
% Full grid for VFS solver format
% Selected grid slices for quick visualizaiton in tecplot

diary grid_matlab_output.txt
clear all
double2 = 0; % double the grid points
coneCylinderFlag = 1; % Cylinder is attached at the base of the cone
onlyCylinderFlag = 0; % only cylinder - no cone
smoothingFlag = 1; % smoothing is applied for better grid continuity
moreThanNintyFlag = 0; % if angle of attack is more than 90 degree (not tested)
AoA = 45; % Angle of attack
rotation = false; % Rotation on x-axis
writeTecPlotFlag = 0; % Write in tecplot format
writeGridDotDatFlag = 1; % write the huge grid.dat file
dispSmoothGrid = 0; % Plot the effect of smooth function

Xmin = -3.3; % Domain limit
Xmax = -Xmin;
g = 0.5;
if(0)
    X1 = incremental_dx_dy2(Xmin,-g,90,1,0.005);
    X2 = incremental_dx_dy2(-g+0.008,-g+g*2/5,30,0,0.003);
    X3 = incremental_dx_dy2(-g+g*2/5+0.003,0-0.008,26,-1,0.008);
    X4 = incremental_dx_dy2(0+0.008,g-g*2/5-0.008,26,1,0.008);
    X5 = incremental_dx_dy2(g-g*2/5,g-0.003,30,0,0.003);
    X6 = incremental_dx_dy2(g,Xmax,90,-1,0.005);
    XTemp = [X1; X2; X3; X4; X5; X6];
    YTemp = [X1; X2; X3; X4; X5; X6];
else
    nnx1 = 150;
    X1 = incremental_dx_dy2(Xmin,-g,nnx1,1,0.005);
    dx = X1(end)-X1(end-1);
    nnx2 = 50;
    X2 = incremental_dx_dy2(-g+dx,0,nnx2,-1,dx);
    X3 = -X2(end-1:-1:1);
    X4 = -X1(end:-1:1);
    XTemp = [X1; X2; X3; X4;];
    YTemp = [X1; X2; X3; X4;];
end

nx = length(XTemp);
ny = length(YTemp);

X = zeros(nx,ny);
Y = zeros(nx,ny);

for i=1:nx
    for j=1:ny
        X(i,j) = XTemp(i);
        Y(i,j) = YTemp(j);
    end
end

if(double2)
    nPx2 = 2*length(X(:,1))-1;
    nPy2 = 2*length(Y(1,:))-1;
    
    x2Temp = zeros(nPx2,nPy2);
    y2Temp = zeros(nPx2,nPy2);
    
    nPxTotal = length(X(:,1));
    nPyTotal = length(Y(1,:));
    
    for i = 1:nPxTotal-1
        for j = 1:nPyTotal-1
            p1x = X(i,j);
            p2x = X(i+1,j);
            p3x = X(i+1,j+1);
            p4x = X(i,j+1);
            p1y = Y(i,j);
            p2y = Y(i+1,j);
            p3y = Y(i+1,j+1);
            p4y = Y(i,j+1);
            
            ii = i*2-1;
            jj = j*2-1;
            
            x2Temp(ii,jj) = p1x;
            x2Temp(ii+1,jj) = 0.5*(p1x+p2x);
            x2Temp(ii+2,jj) = p2x;
            
            x2Temp(ii,jj+1) = 0.5*(p1x+p4x);
            x2Temp(ii+1,jj+1) = 0.5*(0.5*(p1x+p4x) + 0.5*(p2x+p3x));
            x2Temp(ii+2,jj+1) = 0.5*(p2x+p3x);
            
            x2Temp(ii,jj+2) = p4x;
            x2Temp(ii+1,jj+2) = 0.5*(p3x+p4x);
            x2Temp(ii+2,jj+2) = p3x;
            
            y2Temp(ii,jj) = p1y;
            y2Temp(ii+1,jj) = 0.5*(p1y+p2y);
            y2Temp(ii+2,jj) = p2y;
            
            y2Temp(ii,jj+1) = 0.5*(p1y+p4y);
            y2Temp(ii+1,jj+1) = 0.5*(0.5*(p1y+p4y) + 0.5*(p2y+p3y));
            y2Temp(ii+2,jj+1) = 0.5*(p2y+p3y);
            
            y2Temp(ii,jj+2) = p4y;
            y2Temp(ii+1,jj+2) = 0.5*(p3y+p4y);
            y2Temp(ii+2,jj+2) = p3y;
        end
    end
    clear X Y;
    X = x2Temp;
    Y = y2Temp;
end

nPxTotal = length(X(:,1));
nPyTotal = length(Y(1,:));

figure()
surf(X, Y, zeros(nPxTotal,nPyTotal)); view(2);
title("Grid")
drawCircle(0,0,g);
axis equal
axis([-1 1 -1 1])

%% In the section, some necessary variables are calculated for
% grid deformation in z direction.
% Cone base is at 0 in z-axis
% Tip is at some negative value in z-axis (depends on cone angle)
clear factorAtEachLevel factorAtEachLevelFull zFull z xy xyNew C2
baseRadius = 0.5;
zStart = -1.75; zBase = 0; nPz = 300;
dz = 3.000e-03; % Temporary dz
coneAngle = 11.8425; % cone angle
tipz = -baseRadius/tand(coneAngle);
zend_nearTip = tipz - zStart;
eRadius = abs(zend_nearTip)*tand(coneAngle);
%coneAngle = atand((baseRadius-eRadius)/(zBase-zStart));
factor = baseRadius/eRadius;
gap0e = 2*dz/factor;
z = incremental_dx_dy2(zStart,zBase,nPz,0,gap0e*3); % 0:uniform grid generation
% z = linspace(zStart,zBase,nPz)';
radiusatZSlice = baseRadius + (baseRadius - eRadius) * (z - zBase)/(zBase - zStart);
factorInZdir = radiusatZSlice(end:-1:1)/baseRadius; % based on radius
nPz_noScaling1 = 40; % uniform zone near the tip
nPz_noScaling2 = 200; % gradual increase from tip to the inflow
nPz_noScaling3 = 10; % unifrom zone near the inflow
nPzTotal = nPz + nPz_noScaling1 + nPz_noScaling2 + nPz_noScaling3;

% Grid compression or expansion inside the cone
% grids are pushed close to the surface
% C2 decides what will happen to grids inside the body
% Description about C2 is give in the readme file
baseC2 = 0; % near the base
inbetweenBaseTipC2 = 0;
tipC2 = -0.05; % towards the apex
C2(1:nPz) = linspace(baseC2,inbetweenBaseTipC2,nPz);

factoratReverse2 = flip(factorInZdir(end-nPz_noScaling2+1:end))';
C2Reverse = flip(C2(end-nPz_noScaling2+1:end))';

% uniform layer - somewhere inside the body upto the tip
if(nPz_noScaling1 ~= 0)
    gap1 = 2*dz/factor;
    gap1s = z(2)-z(1);
    end1 = gap1s*(nPz_noScaling1+1);
    factorInZdir = [factorInZdir' linspace(factorInZdir(end),factorInZdir(end),nPz_noScaling1)];
    C2 = [C2 linspace(C2(end),0.5*tipC2,nPz_noScaling1)];
    z = [linspace(z(1)-end1, z(1)-gap1s, nPz_noScaling1) z'];
    gap1e = z(end-nPz_noScaling1+2)-z(end-nPz_noScaling1+1);
end

% Gradual expansion layer towards the inflow
if(nPz_noScaling2 ~= 0)
    %reverse 1st layer
    gap2 = 2*dz/factor;
    gap2s = z(2)-z(2-1);
    end2 = 1.25*gap2s*(nPz_noScaling2+1);
    %factorInZdir = [factorInZdir linspace(factorInZdir(end),factorInZdir(end),nPz_noScaling2)];
    factorInZdir = [factorInZdir factoratReverse2];
    C2 = [C2  linspace(C2(end),2.5*tipC2,nPz_noScaling2)];
    %    z = [z linspace(z(end) + gap2, z(end) + end2, nPz_noScaling2)];
    newZ2 = incremental_dx_dy2(z(1) - end2, z(1) - gap2s,nPz_noScaling2,1,gap2s); z = [newZ2' z];
    gap2e = z(end-nPz_noScaling2+2)-z(end-nPz_noScaling2+1);
end

% Uniform zone near the inflow
if(nPz_noScaling3 ~= 0)
    gap3 = 4*dz/factor;
    %gap3 = gap2e;
    gap3s = z(2)-z(1);
    %end3 = 10*dz/factor*(nPz_noScaling3+1);
    factorInZdir = [factorInZdir linspace(factorInZdir(end),factorInZdir(end),nPz_noScaling3)];
    C2 = [C2 linspace(C2(end),C2(end),nPz_noScaling3)];
    %factorInZdir = [factorInZdir factoratReverse3];
    %z = [z linspace(z(end) + gap3, z(end) + end3, nPz_noScaling3)];
    end3 = 1.05*gap3s*(nPz_noScaling3+1); newZ = incremental_dx_dy2(z(1) - end3, z(1) - gap3s,nPz_noScaling3,1,gap3s); z = [newZ' z];
    gap3e = z(end-nPz_noScaling3+2)-z(end-nPz_noScaling3+1);
end

factorInZdir(1:nPzTotal) = factorInZdir(nPzTotal:-1:1);
factorInZdir(1:nPzTotal-nPz+nPz_noScaling3) = smooth(factorInZdir(1:nPzTotal-nPz+nPz_noScaling3),nPz_noScaling3);
C2 = movmean(C2,nPz_noScaling3-1);
C2(1:nPzTotal) = C2(nPzTotal:-1:1);


figure;
plot(z);
hold on;
title("grid points in z direction")
figure;
hold on
plot(z,factorInZdir);
plot(z,C2);
legend('factor','C2');
saveas(gcf,("values_"+mfilename),'png')

%% Transformation is applied on the Cartesian grid to conform to the surface
clear r1
rmax = 3.1; % upto this radius on xy plane compression will take place
% r1 decides upto what radius - which point is considered inside the cone
% and which point is considered outside the cone
r1(1:nPz) = 0.5;
r1(nPz+1 : nPzTotal) = linspace(0.5, 0.51, nPzTotal - nPz);

C2(1:nPzTotal) = C2(nPzTotal:-1:1); %??

sCo = 0.001;
eCo = 0.8;
beta = 0;
xyz = zeros(3,nPxTotal,nPyTotal,nPzTotal);
%z(1:nPzTotal) = z(nPzTotal) - z(nPzTotal:-1:1);

% Controls the size of the grid very very close to the surface
kmax = 1; % corresponding value near the rmax
kmin = 0.95; % corresponding value near the surface (rmin)

% use max function inside the loops to find this value, maximum distance 
% of a point that is inside rmax 
RmaxValue = 3.08; 

% a parabolic function that compress or expand the points in r1 and rmax
syms p q r l m n
A1 = [p^2 p 1;
    q^2 q 1;
    r^2 r 1];
B1 = [l; m; n];
X1 = linsolve(A1,B1);
m1 = 0.9;
m2 = 1.3;
cPositionOfPeak = 0.5;
q1 = r1(1) + cPositionOfPeak * (rmax - r1(1)); % location of the peak of that parabola
xParabola = double(subs(X1,{p,q,r,l,m,n},{r1(1),q1,rmax,m2,m1,1}));
for k = 1:nPzTotal
    for j = 1:nPyTotal
        for i = 1:nPxTotal
            rtemp = sqrt(X(i,j)^2 + Y(i,j)^2);
            if(rtemp <= r1(k)) % inside the cone
                newr = rtemp * exp(C2(k)*(r1(k)-rtemp)/r1(k));
                angle = atan2(Y(i,j),X(i,j));
                myX = newr*cos(angle);
                myY = newr*sin(angle);
                xyz(1,i,j,k) = myX*factorInZdir(k);
                xyz(2,i,j,k) = myY*factorInZdir(k);
            elseif(rtemp > r1(k) && rtemp <= rmax)
                % the idea is to adjust the radius of all the grid points
                % while keeping the angle same for each point 
                % on the Cartesian coordinate. 
                angle = atan2(Y(i,j),X(i,j));
                dist2 = rmax - r1(k);
                rmin = factorInZdir(k)*r1(k);
                dist1 = rmax - rmin;
                %rtemp = r1(k) + (rtemp - r1(k)) * (2 - exp(C2(k)*((rtemp-r1(k))/(rmax - r1(k)))^1));
                %rtemp = r1(k) + (rtemp - r1(k)) * (kmin + (kmax-kmin)*(rtemp-r1(k))/(rmax-r1(k)));
                yV = xParabola(1)*rtemp^2 + xParabola(2)*rtemp + xParabola(3);
                rtemp = r1(k) + (rtemp - r1(k)) * yV;
                dr = rtemp - r1(k);
                newDr = dr * dist1/dist2;
                radiusLocal = factorInZdir(k)*r1(k);
                kvalue = kmin + (kmax-kmin) * (rtemp - r1(k)) / (RmaxValue - r1(k));
                %kvalue = 1;
                %RmaxValue = max(rtemp, RmaxValue);
                myCo = 1 + (sCo - eCo) * (radiusLocal - kvalue*eRadius)/(r1(k) - kvalue*eRadius);
                newDr = newDr - myCo*newDr*(1-newDr/dist1)^1;
                newr = rmin + newDr;
                xyz(1,i,j,k) = newr*cos(angle);
                xyz(2,i,j,k) = newr*sin(angle);
            elseif(rtemp >= rmax) % no transformation
                xyz(1,i,j,k) = X(i,j);
                xyz(2,i,j,k) = Y(i,j);
            end
            xyz(3,i,j,k) = z(k);
        end
    end
    if(0)
        figure()
        xplot(:,:) = xyz(1,:,:,k);
        yplot(:,:) = xyz(2,:,:,k);
        surf(xplot, yplot, zeros(nPxTotal,nPyTotal)+xyz(3,1,1,k)); view(2);
        title("Grid")
        %drawCircle(0,0,0.5);
        axis equal
        %axis([-1 1 -1 1])
    end
end
%% Will add grid points that conform to a cylinder attached at the base of cone
if(coneCylinderFlag) %Cone Cylinder Config
    nPzTotalCylinder = 43;
    dzConeCyl = z(end)-z(end-1);
    lengthOfCylinder = dzConeCyl*nPzTotalCylinder;
    zCyl = linspace(z(end)+lengthOfCylinder/nPzTotalCylinder,z(end)+lengthOfCylinder,nPzTotalCylinder+1);
    xyTemp = xyz;
    %clear xyz;
    
    for k = 1:nPzTotalCylinder
        for i = 1:nPxTotal
            for j = 1:nPyTotal
                xyz(1,i,j,k+nPzTotal) = xyz(1,i,j,nPzTotal);
                xyz(2,i,j,k+nPzTotal) = xyz(2,i,j,nPzTotal);
                xyz(3,i,j,k+nPzTotal) = zCyl(k);
            end
        end
    end
    nPzTotal = length(xyz(1,1,1,:));
end
%% If there is no cone, only cylinder
if (onlyCylinderFlag)
    iSlice = 440; % at what z-location of the cone 
    cylinderradius = factorInZdir(iSlice)*baseRadius;
    disp(['Creating grid for Cylinder of radius: ',num2str(cylinderradius)])
    disp("Achtung: writing grid.dat for cylinder only")
    xyTemp = xyz; clear xyz;
    nPzForCylinderOnly = 20;
    cylinderStart = -0.001;
    lengthOfCylinder = 0.18*cylinderradius;
    zCyl = linspace(cylinderStart,cylinderStart+lengthOfCylinder,nPzForCylinderOnly);
    for k = 1:nPzForCylinderOnly
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                xyz(1,i,j,k) = xyTemp(1,i,j,iSlice);
                xyz(2,i,j,k) = xyTemp(2,i,j,iSlice);
                xyz(3,i,j,k) = zCyl(k);
            end
        end
    end
    nPxTotal = length(xyz(1,:,1,1));
    nPyTotal = length(xyz(1,1,:,1));
    nPzTotal = length(xyz(1,1,1,:));
end

%% Smoothing the grid points - Another better method is applied
% in the next section
tic
if(moreThanNintyFlag==1 && AoA >130)
    nGridSmoothing = 100;
    nGridSmoothing2 = 100;
    nGridSmoothing3 = 30;
else
    nGridSmoothing = 50;
    nGridSmoothing2 = 10;
end

xyNonSmooth = xyz;

if(0) % the later one is more intuitive
    %smoothing on x values
    for k = 1:nPzTotal
        for i = 1:nPxTotal
            xyz(1,i,:,k) = smooth(xyNonSmooth(1,i,:,k),nGridSmoothing);
        end
    end
    for k = 1:nPzTotal
        for i = 1:nPxTotal
            xyz(2,i,:,k) = smooth(xyNonSmooth(2,i,:,k),nGridSmoothing);
        end
    end
    %smoothing on y values
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            xyz(1,:,j,k) = smooth(xyNonSmooth(1,:,j,k),nGridSmoothing);
        end
    end
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            xyz(2,:,j,k) = smooth(xyNonSmooth(2,:,j,k),nGridSmoothing);
        end
    end
    
    %smoothing on z values
    for i = 1:nPxTotal
        for j = 1:nPyTotal
            xyz(3,i,j,:) = smooth(xyNonSmooth(3,i,j,:),nGridSmoothing);
        end
    end
    
    % 2nd smoothing
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            xyz(2,:,j,k) = smooth(xyNonSmooth(2,:,j,k),nGridSmoothing2);
        end
    end
    for k = 1:nPzTotal
        for i = 1:nPxTotal
            xyz(1,i,:,k) = smooth(xyNonSmooth(1,i,:,k),nGridSmoothing2);
        end
    end
    
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            xyz(1,:,j,k) = smooth(xyNonSmooth(1,:,j,k),nGridSmoothing2);
        end
    end
    for k = 1:nPzTotal
        for i = 1:nPxTotal
            xyz(2,i,:,k) = smooth(xyNonSmooth(2,i,:,k),nGridSmoothing2);
        end
    end
    
    %smoothing on z values
    for i = 1:nPxTotal
        for j = 1:nPyTotal
            xyz(3,i,j,:) = smooth(xyNonSmooth(3,i,j,:),nGridSmoothing2);
        end
    end
    
    for k = 1:nPzTotal
        for i=1:nPxTotal
            xyz(1,i,:,k) = smooth(xyz(1,i,:,k),50);
            yLine(:) = xyz(2,i,:,k);
            %plot((xLine),(yLine),'r')
            plot(smooth(xLine,50),smooth(yLine,50),'r')
        end
    end
    
    for k = 1:nPzTotal
        for j=1:nPyTotal
            xLine(:) = xyz(1,:,j,k);
            yLine(:) = xyz(2,:,j,k);
            %plot((xLine),(yLine),'b')
            plot(smooth(xLine,50),smooth(yLine,50),'b')
        end
    end
    
    if(moreThanNintyFlag==1 && AoA >130)
        % 2nd smoothing
        for k = 1:nPzTotal
            for j = 1:nPyTotal
                xyz(2,:,j,k) = smooth(xyNonSmooth(2,:,j,k),nGridSmoothing3);
            end
        end
        for k = 1:nPzTotal
            for i = 1:nPxTotal
                xyz(1,i,:,k) = smooth(xyNonSmooth(1,i,:,k),nGridSmoothing3);
            end
        end
        
        for k = 1:nPzTotal
            for j = 1:nPyTotal
                xyz(1,:,j,k) = smooth(xyNonSmooth(1,:,j,k),nGridSmoothing3);
            end
        end
        for k = 1:nPzTotal
            for i = 1:nPxTotal
                xyz(2,i,:,k) = smooth(xyNonSmooth(2,i,:,k),nGridSmoothing3);
            end
        end
    end
elseif(0) % find a better one below
    xyz = xyNonSmooth;
    %  On xz plane
    for j=1:nPyTotal
        for k=1:nPzTotal
            xyz(1,:,j,k) = smooth(xyNonSmooth(1,:,j,k),20);
        end
        for i=1:nPxTotal
            xyz(1,i,j,1:nPzTotal) = smooth(xyz(1,i,j,1:nPzTotal),50);
        end
    end
    
    % On yz
    for i=1:nPxTotal
        for k=1:nPzTotal
            xyz(2,i,:,k) = smooth(xyz(2,i,:,k),20);
        end
        for j=1:nPyTotal
            xyz(2,i,j,1:nPzTotal) = smooth(xyz(2,i,j,1:nPzTotal),50);
        end
    end
    
    % On xy plane
    for k=1:nPzTotal
        for i=1:nPxTotal
            xyz(1,i,:,k) = smooth(xyz(1,i,:,k),50);
            xyz(2,i,:,k) = smooth(xyz(2,i,:,k),50);
        end
        for j=1:nPyTotal
            xyz(1,:,j,k) = smooth(xyz(1,:,j,k),50);
            xyz(2,:,j,k) = smooth(xyz(2,:,j,k),50);
        end
    end
end % if(smoothingFlag==1)
toc % checking how long the smoothing takes

if(onlyCylinderFlag==1) % swaping coordinates x and z
    xyNew = zeros(3,nPzTotal,nPyTotal,nPxTotal);
    for i = 1:nPxTotal
        for j = 1:nPyTotal
            for k = 1:nPzTotal
                xyNew(2,k,j,i) = xyz(2,i,j,k); % nothing on y
            end
        end
    end
    
    for i = 1:nPxTotal
        for j = 1:nPyTotal
            for k = 1:nPzTotal
                xyNew(3,k,j,i) = xyz(1,i,j,k); % new z are found from old x
            end
        end
    end
    
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                xyNew(1,k,j,i) = xyz(3,i,j,k);  % new x are found from old z
            end
        end
    end
    
    clear xyz;
    xyz = xyNew;
    nPxTotal = length(xyz(1,:,1,1));
    nPyTotal = length(xyz(1,1,:,1));
    nPzTotal = length(xyz(1,1,1,:));
end % onlyCylinderFlag


if(moreThanNintyFlag==1) % swaping coordinates
    xyNew = zeros(3,nPxTotal,nPyTotal,nPzTotal);
    
    for k = 1:nPzTotal
        xyNew(1,:,:,nPzTotal-k+1) = xyz(1,:,:,k);
    end
    for k = 1:nPzTotal
        xyNew(2,:,:,nPzTotal-k+1) = xyz(2,:,:,k);
    end
    for k = 1:nPzTotal
        xyNew(3,:,:,nPzTotal-k+1) = -1*xyz(3,:,:,k);
    end
    
    clear xyz;
    xyz = xyNew;
    nPxTotal = length(xyz(1,:,1,1));
    nPyTotal = length(xyz(1,1,:,1));
    nPzTotal = length(xyz(1,1,1,:));
end

%%  Best way to smooth the grid points : check how smoothing effects the grid
% dispSmoothGrid is for debugging purpose
% On xz plane
if(smoothingFlag==1)
    xyz = xyNonSmooth;
    for j=1:nPyTotal
        for k=1:nPzTotal
            xyz(1,:,j,k) = smooth(xyNonSmooth(1,:,j,k),20);
        end
        for i=1:nPxTotal
            xyz(1,i,j,1:nPzTotal) = smooth(xyz(1,i,j,1:nPzTotal),50);
        end
    end
    if(dispSmoothGrid)
        figure; hold on;
        clear xLine yLine zLine
        j = floor(nPyTotal/2);
        for i=1:nPxTotal
            zLine(:) = xyz(3,i,j,1:nPzTotal);
            xLine(:) = xyz(1,i,j,1:nPzTotal);
            plot((xLine),zLine,'r')
        end
    end
    
    % On yz
    for i = 1:nPxTotal
        for k=1:nPzTotal
            xyz(2,i,:,k) = smooth(xyz(2,i,:,k),20);
        end
        for j=1:nPyTotal
            xyz(2,i,j,1:nPzTotal) = smooth(xyz(2,i,j,1:nPzTotal),50);
        end
    end
    if(dispSmoothGrid)
        figure; hold on;
        clear xLine yLine zLine
        i = floor(nPxTotal/2);
        for j=1:nPyTotal
            zLine(:) = xyz(3,i,j,1:nPzTotal);
            yLine(:) = xyz(2,i,j,1:nPzTotal);
            plot((yLine),(zLine),'r')
        end
    end
    
    if(1)
        tempxy = xyz;
        % On xy plane
        for k = 1:nPzTotal
            for i=1:nPxTotal
                tempxy(2,i,:,k) = smooth(xyz(2,i,:,k),20);
                tempxy(1,i,:,k) = smooth(xyz(1,i,:,k),20);
            end
            for j=1:nPyTotal
                tempxy(1,:,j,k) = smooth(xyz(1,:,j,k),20);
                tempxy(2,:,j,k) = smooth(xyz(2,:,j,k),20);
            end
        end
        xyz = tempxy;
        if(dispSmoothGrid)
            figure; hold on;
            clear xLine yLine zLine
            k = 200;
            for i=1:nPxTotal
                xLine(:) = (xyz(1,i,:,k));
                yLine(:) = (xyz(2,i,:,k));
                plot((xLine),(yLine),'r')
            end
            clear xLine yLine zLine
            for j=1:nPyTotal
                xLine(:) = (xyz(1,:,j,k));
                yLine(:) = (xyz(2,:,j,k));
                plot((xLine),(yLine),'b')
            end
            axis equal
        end
    else
        % On xy plane % looks ok 50 grid for smoothing - but have problems
        figure(13); hold on;
        clear xLine yLine zLine
        k = 200;
        for i=1:nPxTotal
            xLine(:) = smooth(xyz(1,i,:,k),50);
            yLine(:) = smooth(xyz(2,i,:,k),50);
            plot((xLine),(yLine),'r')
        end
        axis equal
        clear xLine yLine zLine
        for j=1:nPyTotal
            xLine(:) = smooth(xyz(1,:,j,k),50);
            yLine(:) = smooth(xyz(2,:,j,k),50);
            plot((xLine),(yLine),'b')
            %plot(smooth(xLine,50),smooth(yLine,50),'b')
        end
        axis equal
    end
    %axis([-0.1 0.1 -2.4 -2.1])
    set(gcf,'Position', [10 10 900 600])
end

%% Checking jacobian
jacobian(1:nPxTotal-1,1:nPyTotal-1,1:nPzTotal-1) = (xyz(1,2:nPxTotal,2:nPyTotal,2:nPzTotal) - xyz(1,1:nPxTotal-1,1:nPyTotal-1,1:nPzTotal-1)) .* ...
    (xyz(2,2:nPxTotal,2:nPyTotal,2:nPzTotal) - xyz(2,1:nPxTotal-1,1:nPyTotal-1,1:nPzTotal-1)) .* ...
    (xyz(3,2:nPxTotal,2:nPyTotal,2:nPzTotal) - xyz(3,1:nPxTotal-1,1:nPyTotal-1,1:nPzTotal-1)) ;
if(any(jacobian(:,:,:)<=0))
    disp("Failed: jacobian negative or zero");
    return;
    beep; pause(1); beep; pause(1); beep; pause(1); beep;
else
    disp("Success: all jacobian are positive");
end
%% 3D representation of the grid
if(1)
    figure(10); hold on;
    for k = 1:5:nPzTotal
        for i=1:5: nPxTotal
            plot3(squeeze(xyz(1,i,:,k)),squeeze(xyz(2,i,:,k)),squeeze(xyz(3,i,:,k)),'r')
        end
    end
    
    for j=1:5: nPyTotal
        for i=1:5: nPxTotal
            plot3(squeeze(xyz(1,i,j,:)),squeeze(xyz(2,i,j,:)),squeeze(xyz(3,i,j,:)),'g')
        end
    end
    
    for k = 1:5:nPzTotal
        for j=1:5: nPyTotal
            plot3(squeeze(xyz(1,:,j,k)),squeeze(xyz(2,:,j,k)),squeeze(xyz(3,:,j,k)),'b')
        end
    end
    axis equal
end

%% 2D visulization of the grid on different orthogonal planes
if(1)
    %On xz plane
    figure(11); hold on;
    clear xLine yLine zLine
    j = floor(nPyTotal/2)+0;
    for i=1:nPxTotal
        zLine(:) = xyz(3,i,j,1:nPzTotal);
        xLine(:) = xyz(1,i,j,1:nPzTotal);
        plot((xLine),zLine,'r')
    end
    axis equal
    title('xz plane at y=mid')
    
    % On yz
    figure(12); hold on;
    clear xLine yLine zLine
    i = floor(nPxTotal/2);
    for j=1:nPyTotal
        zLine(:) = xyz(3,i,j,1:nPzTotal);
        yLine(:) = xyz(2,i,j,1:nPzTotal);
        plot((yLine),(zLine),'r')
    end
    axis equal
    title('yz plane at x=mid')
    
    
    % On xy plane
    figure(13); hold on;
    clear xLine yLine zLine
    k = 200;
    for i=1:nPxTotal
        xLine(:) = (xyz(1,i,:,k));
        yLine(:) = (xyz(2,i,:,k));
        plot((xLine),(yLine),'r')
    end
    axis equal
    clear xLine yLine zLine
    for j=1:nPyTotal
        xLine(:) = (xyz(1,:,j,k));
        yLine(:) = (xyz(2,:,j,k));
        plot((xLine),(yLine),'b')
    end
    axis equal
    title("xy plane at z="+xyz(3,1,1,k))
    
    
    % On xz plane
    figure(14); hold on;
    clear xLine yLine zLine
    j = 100;
    for k=1:nPzTotal
        xLine(:) = (xyz(1,:,j,k));
        zLine(:) = (xyz(3,:,j,k));
        plot((xLine),(zLine),'r')
    end
    axis equal
    title("xz plane at y="+xyz(2,1,j,1))
    
    %axis([-0.1 0.1 -2.4 -2.1])
    %set(gcf,'Position', [10 10 900 600])
end

%% Writing down some info about the grid
datafiles = dir('yy_TecPlot3D*.dat');
lastFilenumber = length(datafiles);
if (lastFilenumber ~= 0)
    ID = sscanf(datafiles(lastFilenumber).name, 'yy_TecPlot3D.%d') + 1;
else
    ID = 1;
end
IDstr = sprintf('%07d',ID);
fname = ['yy_Grid_info.',IDstr,'.txt'];
fid = fopen(fname,'W');
fprintf(fid,'%s,\n',datetime);
fprintf(fid,'baseRadius = %f, eRadius = %f \n',baseRadius, eRadius);
fprintf(fid,'nPxTotal = %d, nPyTotal = %d, nPzTotal = %d \n', nPxTotal, nPyTotal, nPzTotal);
fprintf(fid,'zStart = %f; zBase = %f; nPz = %d; dz = %f; \n', zStart, zBase, nPz, dz);
fprintf(fid,'nPz_noScaling1 = %d; nPz_noScaling2 = %d; nPz_noScaling3 = %d \n', nPz_noScaling1, nPz_noScaling2, nPz_noScaling3);
fprintf(fid,'gap0e = %f, gap1s = %f, gap1 = %f, gap1e = %f \n', gap0e, gap1s, gap1, gap1e);
if(nPz_noScaling2 ~= 0)
    fprintf(fid,'gap2s = %f, gap2 = %f, gap2e = %f \n', gap2s, gap2, gap2e);
end
if(nPz_noScaling3 ~= 0)
    fprintf(fid,'gap3s = %f, gap3 = %f, gap3e = %f \n', gap3s, gap3, gap3e);
end
fprintf(fid,'rmax = %f; r1(1) = %f; sCo = %f, eCo = %f , beta = %f \n', rmax, r1(1), sCo, eCo, beta);
fprintf(fid,'kmax = %f; kmin = %f; RmaxValue = %f, \n m1 = %f , m2 = %f, cPositionOfPeak = %f \n', kmax, kmin, RmaxValue, m1, m2, cPositionOfPeak);

if(onlyCylinderFlag)
    fprintf(fid,'onlyCylinderFlag = %d, nPzForCylinderOnly = %d, cylinderStart = %f lengthOfCylinder = %f\n', onlyCylinderFlag, nPzForCylinderOnly, cylinderStart, lengthOfCylinder);
end
if(smoothingFlag)
    fprintf(fid,'smoothingFlag = %d, nGridSmoothing = %d, nGridSmoothing2 = %d \n', smoothingFlag, nGridSmoothing, nGridSmoothing2);
end
if(coneCylinderFlag)
    fprintf(fid,'coneCylinderFlag = %d, nPzTotalCylinder = %d, lengthOfCylinder = %f \n', coneCylinderFlag, nPzTotalCylinder, lengthOfCylinder);
end
fclose(fid);
disp('done writing in Tachpolot 3d')

%% The following one can be opened in visit - fewer grid points
fname=['yy_TecPlot3D.',IDstr,'.dat'];
disp(['writing: ', fname]);
fid = fopen(fname,'W');
fprintf(fid,'TITLE = "Cylindrical Flow Velocity Data"\n');
fprintf(fid,'Variables= "x","y","z","u"\n');

clear slicesZ
slicesZ = [1 nPz_noScaling3 nPz_noScaling3+nPz_noScaling2 nPz_noScaling3+nPz_noScaling2+nPz_noScaling1 floor((nPz_noScaling3+nPz_noScaling2+nPz_noScaling1+nPzTotal)/2) nPzTotal];
if(coneCylinderFlag)
    slicesZ = [1 nPz_noScaling3 nPz_noScaling3+nPz_noScaling2 nPz_noScaling3+nPz_noScaling2+nPz_noScaling1 floor((nPz_noScaling3+nPz_noScaling2+nPz_noScaling1+nPzTotal)/2) nPzTotal-nPzTotalCylinder nPzTotal];
end

if(onlyCylinderFlag == 1)
    nPzTotalTec = length(xyz(1,1,1,:));
elseif(moreThanNintyFlag == 1)
    slicesZ = nPzTotal - slicesZ + 1;
    nPzTotalTec = length(slicesZ);
else
    nPzTotalTec = length(slicesZ);
end

fprintf(fid,'ZONE I=%d, J=%d, K=%d, F=POINT \n',nPxTotal, nPyTotal, nPzTotalTec);

if(onlyCylinderFlag == 1)
    for k = 1:nPzTotalTec
        for j = 1:nPyTotal
            for i = 1 : nPxTotal
                fprintf(fid,'%15e %15e %15e 1\n',xyz(1,i,j,k),xyz(2,i,j,k),xyz(3,i,j,k));
            end
        end
    end
elseif(moreThanNintyFlag == 1)
    for l = 1:nPzTotalTec
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                k = slicesZ(l);
                fprintf(fid,'%15e %15e %15e 1\n',xyz(1,i,j,k),xyz(2,i,j,k),xyz(3,i,j,k));
            end
        end
    end
else
    for l = 1:nPzTotalTec
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                k = slicesZ(l);
                fprintf(fid,'%15e %15e %15e 1\n',xyz(1,i,j,k),xyz(2,i,j,k),xyz(3,i,j,k));
            end
        end
    end
end
fclose(fid);
disp(['Done writing: ', fname]);
% beep; pause(0.6); beep; pause(0.6); beep; pause(0.6); beep;
%% VFS solver format
if(writeGridDotDatFlag==1)
    dim = 3;
    nPxTotal = length(xyz(1,:,1,1));
    nPyTotal = length(xyz(1,1,:,1));
    nPzTotal = length(xyz(1,1,1,:));
    fname = sprintf('grid.dat');
    fid = fopen(fname,'W');
    fprintf(fid,'1 \n');
    fprintf(fid,'%d  %d  %d \n', nPxTotal, nPyTotal, nPzTotal);
    
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                fprintf(fid,'%15e ',xyz(1,i,j,k));
            end
        end
    end
    fprintf(fid,'\n');
    
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                fprintf(fid,'%15e ',xyz(2,i,j,k));
            end
        end
    end
    fprintf(fid,'\n');
    
    for k = 1:nPzTotal
        for j = 1:nPyTotal
            for i = 1:nPxTotal
                fprintf(fid,'%15e ',xyz(3,i,j,k));
            end
        end
    end
    
    fclose(fid);
end
%% Writing down necessary info about the grid
if(0)    
    IDstr = sprintf('%07d',0);
    fname = ['Grid_info.',IDstr,'.txt'];
    fid = fopen(fname,'W');
    fprintf(fid,'%s,\n',datetime);
    fprintf(fid,'baseRadius = %f, eRadius = %f \n',baseRadius, eRadius);
    fprintf(fid,'nPxTotal = %d, nPyTotal = %d, nPzTotal = %d \n', nPxTotal, nPyTotal, nPzTotal);
    fprintf(fid,'zStart = %f; zBase = %f; nPz = %d; dz = %f; \n', zStart, zBase, nPz, dz);
    fprintf(fid,'nPz_noScaling1 = %d; nPz_noScaling2 = %d; nPz_noScaling3 = %d \n', nPz_noScaling1, nPz_noScaling2, nPz_noScaling3);
    fprintf(fid,'gap0e = %f, gap1s = %f, gap1 = %f, gap1e = %f \n', gap0e, gap1s, gap1, gap1e);
    if(nPz_noScaling2 ~= 0)
        fprintf(fid,'gap2s = %f, gap2 = %f, gap2e = %f \n', gap2s, gap2, gap2e);
    end
    if(nPz_noScaling3 ~= 0)
        fprintf(fid,'gap3s = %f, gap3 = %f, gap3e = %f \n', gap3s, gap3, gap3e);
    end
    fprintf(fid,'rmax = %f; r1(1) = %f; sCo = %f, eCo = %f , beta = %f \n', rmax, r1(1), sCo, eCo, beta);
    fprintf(fid,'kmax = %f; kmin = %f; RmaxValue = %f, \n m1 = %f , m2 = %f, cPositionOfPeak = %f \n', kmax, kmin, RmaxValue, m1, m2, cPositionOfPeak);
    
    if(onlyCylinderFlag)
        fprintf(fid,'onlyCylinderFlag = %d, nPzForCylinderOnly = %d, cylinderStart = %f lengthOfCylinder = %f\n', onlyCylinderFlag, nPzForCylinderOnly, cylinderStart, lengthOfCylinder);
    end
    if(smoothingFlag)
        fprintf(fid,'smoothingFlag = %d, nGridSmoothing = %d, nGridSmoothing2 = %d \n', smoothingFlag, nGridSmoothing, nGridSmoothing2);
    end
    if(coneCylinderFlag)
        fprintf(fid,'coneCylinderFlag = %d, nPzTotalCylinder = %d, lengthOfCylinder = %f \n', coneCylinderFlag, nPzTotalCylinder, lengthOfCylinder);
    end
    
    fclose(fid);
end

%% Rotation
if(0)
    if(rotation)
        angleTheta = 90;
        Rmat = [1           0               0;        ...
            0 cosd(angleTheta)  sind(angleTheta); ...
            0 -sind(angleTheta) cosd(angleTheta)];
        
        for k = 1:nPzTotal
            for i = 1:nPxTotal
                for j = 1:nPyTotal
                    xyNew(:,i,j,k) = Rmat*xyz(:,i,j,k);
                end
            end
        end
        xyz = xyNew;
    end
end

%% write data in tecplot format
if(writeTecPlotFlag)
    dim = 3;
    fname = sprintf('gridTecplot.dat');
    fid = fopen(fname,'wt');
    fprintf(fid,'TITLE = "Flow Velocity Data"\n');
    fprintf(fid,'Variables= "x","y","z","u"\n');
    fprintf(fid,'ZONE I=%d, J=%d, K=%d, F=POINT \n',nPxTotal, nPyTotal, nPzTotal);
    
    for k = 1 : nPzTotal
        for j = 1 : nPyTotal
            for i = 1: nPxTotal
                fprintf(fid,'%15e %15e %15e 1\n',xyz(1,i,j,k),xyz(2,i,j,k),xyz(3,i,j,k));
            end
        end
    end
    fprintf(fid,'\n');
    
    fclose(fid);
end