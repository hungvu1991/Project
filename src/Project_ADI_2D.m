%Hung Vu
%UHID: 1152030
%Project Bc1-4 with ADI Method with error
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

%% Change 
dx=x(2)-x(1); %Change in x-direction
dy=y(2)-y(1); %Change in y-direction
dt=.125*(dx^2+dy^2); %Change in t-direction with Time step
lambda=dt/dx^2; %Setting lambda for better computation;
lambda2=2*lambda; %Setting lambda for better computation;
XX=2*(1-lambda); %Better outside the loop

%% Boundary Conditions
faby=(by-ay)^2*cos(pi*by/ay); %Given for the top boundary condition
gaby=by*(by-ay)^2; %Given for the top boundary condition
Left=(y-ay).^2.*cos(pi*y/ay); %Left Boundary Conditions
Right=y.*(y-ay).^2; %Right Boundary Conditions
Top=(faby+(x-ax)./(bx-ax)*(gaby-faby)); %Top Boundary Conditions
Bottom=zeros(1,length(x)); % Bottom Boundary Condition

%% Variables
L1=length(y)-2; %Setting up Length
L2=L1+1; %Setting up Length
Bas=2:L1+1; %Setting up Length

%% Tridiagonal matrix algorithm setup
dia1=2*(1+lambda)*ones(1,L1); %Creating a row to put in the diagonal
dia2=-lambda*ones(1,L1-1); %Creating a row to put in the diagonal
Mat=diag(dia1); %Adding diagonal into matrix
Mat=diag(dia2,1)+Mat; %Adding diagonal into matrix
Mat=diag(dia2,-1)+Mat; %Adding diagonal into matrix
mat=Mat; %Same matrix
mat(L1,L1+1)=Mat(1,2); %Adding more entries
mat(L1+1,[L1,L1+1])=[-lambda2,Mat(L1,L1)]; %Adding more entries

%% Adding boundary conditions into a U matrix.
U=zeros(Nx-2,Nx-2); % Initial condition with boundary condition
TRU=zeros(Nx,Ny); %Creating matrix with 0's
TRU(L1+2:-1:1,1)=Left; %Adding left boundary condition
TRU(1,1:1:L1+2)=Top; %Adding Top boundary condition
TRU(L1+2:-1:1,Nx)=Right; %Adding left boundary condition

%% Checkpoint
if exist( 'checkpointEX.mat','file' ) % If a checkpoint file exists, load it
    load('checkpointEX.mat')
end

Count=0;
%% For loop
Poi=1;
while Poi < 10000

    %% Creating the first row with boundary conditions
    A(1,1:L1) = lambda*Top(Bas)+XX*U(1,1:L1)+lambda*U(2,1:L1);
    A(1) = A(1) + lambda*Left(2);
    A(L1) = A(L1) + lambda*Right(2);
    U(1,:) = Mat\A';
    
    %% Creating the rest of the rows with boundary conditions
    AA(1:L1-2,1:L1) = lambda*U(1:L1-2,1:L1)+XX*U(2:L1-1,1:L1)+lambda*U(3:L1,1:L1);
    AA(:,1) = AA(:,1) + lambda*Left(3:L1)';
    AA(:,L1) = AA(:,L1) + lambda*Right(3:L1)';
    U(2:L1-1,:) = (Mat\AA')';
    
    %% Creating the last row with boundary conditions
    A(1,1:L1) = lambda*U(L1-1,1:L1)+XX*U(L1,1:L1)+lambda*Bottom(Bas);
    A(1) = A(1) + lambda*Left(L2);
    A(L1) = A(L1) + lambda*Right(L2);
    U(L1,:) = Mat\A';
    
    %% Creating the bottom condtion
    A(1,1:L1) = lambda2*U(L1,1:L1)+XX*Bottom(Bas);
    A(1) = A(1) + lambda*Left(L2);
    A(L1) = A(L1) + lambda*Right(L2);
    Bottom(2:L1+1) = Mat\A';
    
    %% Creating the first column with boundary conditions
    B(1,1:L1) = lambda*Left(Bas)'+XX*U(1:L1,1)+lambda*U(1:L1,2);
    B(1) = B(1) + lambda*Top(2);
    B(L2) = lambda*Left(L1+2)+XX*Bottom(2)+lambda*Bottom(3);
    T_dummy = (mat\B')';
    U(:,1) = T_dummy(1:L1);
    Bottom(2) = T_dummy(L2);
    
    %% Creating the rest of the columns with boundary conditions
    j = 2:L1-1;
    BB(1:L1,1:L1-2) = lambda*U(1:L1,1:L1-2)+XX*U(1:L1,j)+lambda*U(1:L1,3:L1);
    BB(1,:) = BB(1,:) + lambda*Top(3:L1);
    BB(L1+1,:) = lambda*Bottom(2:L1-1)+XX*Bottom(3:L1)+lambda*Bottom(4:L1+1);
    T_dummy = mat\BB;
    U(:,2:L1-1) = T_dummy(1:L1,1:L1-2);
    Bottom(3:L1) = T_dummy(L2,:);
    
    %% Creating the last column with boundary conditions
    B(1,1:L1) = lambda*U(1:L1,L1-1)+XX*U(1:L1,L1)+lambda*Right(Bas)';
    B(1) = B(1) + lambda*Top(L1+1);
    B(L2) = lambda*Bottom(L1)+XX*Bottom(L2)+lambda*Right(L1+2);
    T_dummy = (mat\B')';
    U(:,L1) = T_dummy(1:L1);
    Bottom(L2) = T_dummy(L2);
    
    %% Updating the true values
    TRU(Nx,L1+2:-1:1) = Bottom;
    TRU(2:L1+1,L1+1:-1:2) = U;
    
    if (mod(Poi,10^4) == 0)
        disp('Iteration passed:')
        disp(Poi)
    end
    
    Count=Count+1;  %Increase the count
    if mod(Count,1000)==0   %Save checkpoint file every 1000 iterations
        save('checkpointEX.mat'); %Save the file
    end     %Close if loop
    
    Poi = Poi+1; % Going to the next iterations
end

clc
Time_SS = dt*Poi %Finding the steady state value


%% 5% True value test
clearvars -except TRU time_SS

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
[X,Y]=meshgrid(-x,-y); %Creating a mesh

%% Change
dx=x(2)-x(1); %Change in x-direction
dy=y(2)-y(1); %Change in x-direction
dt=.125*(dx^2+dy^2); %Change in t-direction
lambda=dt/dx^2; %Setting lambda for better computation;
lambda2=2*lambda; %Setting lambda for better computation;
XX=2*(1-lambda); %Better outside the loop

%% Boundary Conditions
faby=(by-ay)^2*cos(pi*by/ay); %Given for the top boundary condition
gaby=by*(by-ay)^2; %Given for the top boundary condition
Left=fliplr((y-ay).^2.*cos(pi*y/ay)); %Left Boundary Conditions
Right=fliplr(y.*(y-ay).^2); %Right Boundary Conditions
Top=(faby+(x-ax)./(bx-ax)*(gaby-faby)); %Top Boundary Conditions
Bottom=zeros(1,length(x)); % right side

%% varaibles used outside while loop to speed up process:
L1=length(y)-2; %Setting up Length
L2=L1+1; %Setting up Length
Bas=2:L1+1; %Setting up Length

%% Tridiagonal matrix algorithm setup
dia1=2*(1+lambda)*ones(1,L1); %Creating a row to put in the diagonal
dia2=-lambda*ones(1,L1-1); %Creating a row to put in the diagonal
Mat=diag(dia1); %Adding diagonal into matrix
Mat=diag(dia2,1)+Mat; %Adding diagonal into matrix
Mat=diag(dia2,-1)+Mat; %Adding diagonal into matrix
mat=Mat; %Same matrix
mat(L1,L1+1)=Mat(1,2); %Adding more entries
mat(L1+1,[L1,L1+1])=[-lambda2,Mat(L1,L1)]; %Adding more entries

%% Creating a U matrix for an approximation value.
APX=zeros(Nx,Ny); %Creating matrix with 0's
APX(L1+2:-1:1,1)=Left; %Adding left boundary conditions
APX(1,1:1:L1+2)=Top; %Adding top boundary conditions
APX(L1+2:-1:1,Nx)=Right; %Adding right boundary conditions
U=zeros(Nx-2,Nx-2); % inital condition

%% Checkpoint
if exist( 'checkpointEX.mat','file' ) % If a checkpoint file exists, load it
    load('checkpointEX.mat')
end
Count=0;

%% While loop
Poi=1;
error=inf;

while error > 5

    %% Creating the first row with boundary conditions
    A(1,1:L1) = lambda*Top(Bas)+XX*U(1,1:L1)+lambda*U(2,1:L1);
    A(1) = A(1) + lambda*Left(2);
    A(L1) = A(L1) + lambda*Right(2);
    U(1,:) = Mat\A';
    
    %% Creating the rest of the rows with boundary conditions
    AA(1:L1-2,1:L1) = lambda*U(1:L1-2,1:L1)+XX*U(2:L1-1,1:L1)+lambda*U(3:L1,1:L1);
    AA(:,1) = AA(:,1) + lambda*Left(3:L1)';
    AA(:,L1) = AA(:,L1) + lambda*Right(3:L1)';
    U(2:L1-1,:) = (Mat\AA')';
    
    %% Creating the last row with boundary conditions
    A(1,1:L1) = lambda*U(L1-1,1:L1)+XX*U(L1,1:L1)+lambda*Bottom(Bas);
    A(1) = A(1) + lambda*Left(L2);
    A(L1) = A(L1) + lambda*Right(L2);
    U(L1,:) = Mat\A';
    
    %% Creating the bottom boundary condition
    A(1,1:L1) = lambda2*U(L1,1:L1)+XX*Bottom(Bas);
    A(1) = A(1) + lambda*Left(L2);
    A(L1) = A(L1) + lambda*Right(L2);
    Bottom(2:L1+1) = Mat\A';

    %% Creating the first column with boundary conditions
    B(1,1:L1) = lambda*Left(Bas)'+XX*U(1:L1,1)+lambda*U(1:L1,2);
    B(1) = B(1) + lambda*Top(2);
    B(L2) = lambda*Left(L1+2)+XX*Bottom(2)+lambda*Bottom(3);
    T_dummy = (mat\B')';
    U(:,1) = T_dummy(1:L1);
    Bottom(2) = T_dummy(L2);
    
    %% Creating the rest of the columns with boundary conditions
    j = 2:L1-1;
    BB(1:L1,1:L1-2) = lambda*U(1:L1,1:L1-2)+XX*U(1:L1,j)+lambda*U(1:L1,3:L1);
    BB(1,:) = BB(1,:) + lambda*Top(3:L1);
    BB(L1+1,:) = lambda*Bottom(2:L1-1)+XX*Bottom(3:L1)+lambda*Bottom(4:L1+1);
    T_dummy = mat\BB;
    U(:,2:L1-1) = T_dummy(1:L1,1:L1-2);
    Bottom(3:L1) = T_dummy(L2,:);
    
    %% Creating the last column with boundary conditions
    B(1,1:L1) = lambda*U(1:L1,L1-1)+XX*U(1:L1,L1)+lambda*Right(Bas)';
    B(1) = B(1) + lambda*Top(L1+1);
    B(L2) = lambda*Bottom(L1)+XX*Bottom(L2)+lambda*Right(L1+2);
    T_dummy = (mat\B')';
    U(:,L1) = T_dummy(1:L1);
    Bottom(L2) = T_dummy(L2);
    
    %% Updating the approximated values 
    APX(Nx,L1+2:-1:1) = Bottom;
    APX(2:L1+1,L1+1:-1:2) = U;
    
    error = abs(sum(sum(TRU-APX)))/Nx^2*100; %calulating the errors
    
    %% Creating a blot
    drawnow %updates
    Time = Poi*dt; %calulating the Time
    if ( mod(Poi,100) == 0)
        contourf(X,Y,APX)%Creating contour plot
        caxis( [ min(min(APX)), max(max(APX)) ] ) %Creating axis values
        colorbar %adding colorbar to indicate values with color since 2D
        title(sprintf('Time (secs) %10.4f seconds',Time)) %Title
        drawnow %updates
    end
    
    Poi = Poi+1; % next iteration
      
    
    Count=Count+1;  %Increase the count
    if mod(Count,1000)==0   %Save checkpoint file every 1000 iterations
        save('checkpointEX.mat'); %Save the file
    end %close if loop
end

%% final plot
drawnow %updates
Time = Poi*dt; %calulating the Time
time_apprx = Time; %calulating the Time
contourf(X,Y,APX) %calulating the Time
caxis( [ min(min(APX)), max(max(APX)) ] ) %Creating axis values
colorbar %adding colorbar to indicate values with color since 2D
title(sprintf('Time (secs) %10.4f seconds',Time)) %Title
drawnow %updates

disp('Approximate time it took:')
disp(time_apprx)
disp('Error:')
disp(error)
delete('checkpointEX.mat');    %Delete checkpoint file once evertything is complete