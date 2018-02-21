%Phiip Johnson
%Plotter for hal results

clc
clear all


RES = 3306; %number of elements
Ns = 9; %Solution points per element

D = 0.001;

%Get all of the x,y points in the mesh and the initial distribution
RHO_IN = csvread('friendly_rho.csv');
Vx_IN = csvread('friendly_ux.csv');
Vy_IN = csvread('friendly_uy.csv');
P_IN = csvread('friendly_p.csv');

for e = 1:RES
    for g = 1:Ns
        j = (e-1)*Ns + g;
        X(j) = RHO_IN(j,1)/D;
        Y(j) = RHO_IN(j,2)/D;
        RHO(j) = RHO_IN(j,4);
        Vx(j) = Vx_IN(j,1);
        Vy(j) = Vy_IN(j,1); 
        P(j) = P_IN(j,1);
        SPEED(j) = sqrt(Vx(j)^2+Vy(j)^2);
    end
end



n=32; %number of subdivisions for plotting purposes
[xgrid,ygrid]=meshgrid(linspace(min(X),max(X),32*n),linspace(min(Y),max(Y),6*n));
[xgZ,ygZ]=meshgrid(linspace(-12,12,12*n),linspace(min(Y),max(Y),3*n));
RHO_grid = griddata(X,Y,RHO,xgrid,ygrid);
RHO_gZ = griddata(X,Y,RHO,xgZ,ygZ);
Vx_grid = griddata(X,Y,Vx,xgrid,ygrid);
Vx_gZ = griddata(X,Y,Vx,xgZ,ygZ);
Vy_grid = griddata(X,Y,Vy,xgrid,ygrid);
Vy_gZ = griddata(X,Y,Vy,xgZ,ygZ);
P_grid = griddata(X,Y,P,xgrid,ygrid);
P_gZ = griddata(X,Y,P,xgZ,ygZ);

SPEED_gZ = griddata(X,Y,SPEED,xgZ,ygZ);

%figure(1)
%surf(xgrid,ygrid,Vx_grid)
%shading flat
%axis tight



%figure(3)
%quiver3(xgrid,ygrid,zgrid, Vx_grid, Vy_grid, Vz_grid,1.0)

% figure(1)
% surface(xgrid,ygrid,RHO_grid)
% shading flat
% colormap hot
% colorbar
% caxis([.305 .31])
% xlabel('x/D')
% ylabel('y/D')
% %axis square
% pbaspect([40 6 1])
% axis([-32 32 0 6])
% view(-0,90)    

figure(2)
surface(xgZ,ygZ,RHO_gZ)
shading flat
colormap hot
colorbar
caxis([.305 .31])
xlabel('x/D')
ylabel('y/D')
%axis square
pbaspect([4 1 1])
axis([-12 12 0 6])
view(-0,90)    

figure(5)
surface(xgZ,ygZ,Vy_gZ)
shading flat
colormap jet
colorbar
xlabel('x/D')
ylabel('y/D')
pbaspect([4 1 1])
axis([-12 12 0 6])
view(-0,90)    

figure(6)
surface(xgZ,ygZ,SPEED_gZ)
shading flat
colormap jet
colorbar
xlabel('x/D')
ylabel('y/D')
pbaspect([4 1 1])
axis([-12 12 0 6])
view(-0,90)    

figure(3)
surface(xgZ,ygZ,Vx_gZ)
shading flat
colormap jet
colorbar
xlabel('x/D')
ylabel('y/D')
pbaspect([4 1 1])
axis([-12 12 0 6])
view(-0,90)    

figure(4)
surface(xgZ,ygZ,P_gZ)
shading flat
colormap hot
colorbar
xlabel('x/D')
ylabel('y/D')
pbaspect([4 1 1])
axis([-12 12 0 6])
view(-0,90)    