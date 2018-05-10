%Hung Vu
%UHID: 1152030
%Project Bc1-4 with Explicit Method
clear all; clc; close all;

%% domain of interest
ax=-pi; %parameter of x
bx=pi;  %parameter of x
ay=-pi; %parameter of y
by=pi;  %parameter of y

%% Nodes
Nx=100; %Number of Nodes for x-direction
Ny=100; %Number of Nodes for y-direction
x=linspace(ax,bx,Nx);   %Creaing a matrix
y=linspace(ay,by,Ny);   %Creaing a matrix
[X,Y]=meshgrid(x,y);    %Creaing a matrix

%% Change in matrix
dx=x(2)-x(1); %Change in x-direction
dy=y(2)-y(1); %Change in y-direction
dt=min([dx,dy])^2/8; %Change in t-direction with time step
lambda=dt/dx^2 %Setting lambda for better computation; Less than 0.25 stable

%% Checkpoint
if exist( 'checkpointEX.mat','file' ) % If a checkpoint file exists, load it
    load('checkpointEX.mat')
end

%% U matrix
U=zeros(Nx,Ny);  %Creating martrix for initial u conditions
time=0; %Initial time condition
Poi=1000; % Total iterations being tested
Count=0; %Count for checkpoint

%% Boundary considtions with diffusion
for k=1:Poi
    U_New=zeros(Nx,Ny); %Setting new u matrix for diffusion
    for j=2:Ny-1
        for i=2:Nx-1
            U_New(i,j)=(lambda*(-4*U(i,j)+U(i,j-1)+U(i-1,j)+U(i+1,j)+U(i,j+1)))+U(i,j); %Explicit
        end
    end
    
    %% Adding boundary conditions
    U_New(:,1)=fliplr((y-ay).^2.*cos(pi*y/ay)); %Left Boundary Conditions
    U_New(:,Ny)=fliplr(y.*(y-ay).^2); %Right Boundary Conditions
    faby=(by-ay)^2*cos(pi*by/ay); %Given for the top boundary condition
    gaby=by*(by-ay)^2; %Given for the top boundary condition
    U_New(1,:)=faby+(x-ax)./(bx-ax)*(gaby-faby);  %Top Boundary Conditions
    
    %% Adding ghost node(s) for Neumann condition
    for i=2:Nx-1
       U_New(Nx,i)=(lambda*(U(i,j-1)+2*U(i-1,j)+U(i,j+1)))+(1-2*lambda-2*lambda)*U(i,j); %Explicit
    end
   
    %% Converting the u matrix
    U=U_New; 
    time=time+dt;  %Change in time
    U_New1=flipud(U_New); %Flipping the matrix for graph
    
    %% Creating the plot
    if (mod(k,100) == 0)
        contourf(X,Y,U_New1) %Creating contour plot
        caxis([min(min(U_New1)) max(max(U_New1))]) %Creating axis values
        colorbar %adding colorbar to indicate values with color since 2D
        title(sprintf('Time (secs) %10.4f',time)) %Title
        drawnow %updates
    end
    
    Count=Count+1;  %Increase the count
    if mod(Count,1000)==0   %Save checkpoint file every 1000 iterations
        save('checkpointEX.mat'); %Save the file
    end     %Close if loop
    
end

delete('checkpointEX.mat');    %Delete checkpoint file once evertything is complete