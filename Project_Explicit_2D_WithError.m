%Hung Vu
%UHID: 1152030
%Project Bc1-4 with Explicit Method with Error
clear all; clc; close all;

%Domain of interest
ax=-pi; %parameter of x
bx=pi;  %parameter of x
ay=-pi; %parameter of y
by=pi;  %parameter of y

Nx=100; %Number of Nodes for x-direction
Ny=100; %Number of Nodes for y-direction
x=linspace(ax,bx,Nx);   %Creaing a matrix
y=linspace(ay,by,Ny);   %Creaing a matrix
[X,Y]=meshgrid(x,y);    %Creaing a matrix

dx=x(2)-x(1); %Change in x-direction
dy=y(2)-y(1); %Change in y-direction
dt=min([dx,dy])^2/8; %Change in t-direction with time step
lambda=dt/dx^2; %Setting lambda for better computation; Less than 0.25 stable

U=zeros(Nx,Nx);  %Setting initial temperature conditions
time=0; %Setting initial time
err=1; %Setting initial error
Poi=1; % total iterations
XX=zeros(Nx,Nx);

while err > 10^-10
    U_new=zeros(Nx,Nx);
    for j=2:Nx-1
        for i=2:Nx-1
            U_new(i,j)=(lambda*(-4*U(i,j)+U(i,j-1)+U(i-1,j)+U(i+1,j)+U(i,j+1)))+U(i,j); %Explicit
        end
    end
    
    % Adding boundary conditions
    U_new(:,1)=fliplr((y-ay).^2.*cos(pi*y/ay)); %Left Boundary Conditions
    U_new(:,Ny)=fliplr(y.*(y-ay).^2); %Right Boundary Conditions
    faby=(by-ay)^2*cos(pi*by/ay); %Given for the top boundary condition
    gaby=by*(by-ay)^2; %Given for the top boundary condition
    U_new(1,:)=faby+(x-ax)./(bx-ax)*(gaby-faby);  %Top Boundary Conditions
    %%
    %Adding ghost node(s) for Neumann condition
    for i=2:Nx-1
        U_new(Nx,i)=(lambda*(U(i,j-1)+2*U(i-1,j)+U(i,j+1)))+(1-2*lambda-2*lambda)*U(i,j); %Explicit
    end
    
    U=U_new; %Converting the u matrix
    time=time+dt; %Change in time
    U_new=flipud(U_new); %Flipping the matrix for graph
    
    %error check:
    a=mean(mean(U_new));
    b=mean(mean(XX));
    err=abs((a-b)/a*100);
    
    if (mod(Poi,1000) == 0)
        disp('Iterations:')
        disp(Poi)
        disp('Error:')
        disp(err)
    end
    Poi=1+Poi;
    XX=U_new;
    
end