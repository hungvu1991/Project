clear, clc

% MECE 5397
% Bc2-1 Diffusion Equation
% Using ADI Scheme

%% Discretization and Inputing Given Data:

dis = 100; % total number of nodes, dis x dis


x = linspace(-pi,pi,dis);
y = linspace(-pi,pi,dis);
del_x = x(2)-x(1);
del_y = y(2)-y(1);
[X,Y]= meshgrid(x,y);

a_x = -pi; a_y = -pi;
b_x = pi;  b_y = pi;

dt= min([del_x,del_y])^2/8;       %Computes stable time step
lambda_x=dt/del_x^2;
lambda_y=dt/del_y^2;
T=zeros(dis,dis);  %Set initial temperature conditions
time=0;

%%
err = 1;
iter = 1; % total iterations
T_previous = zeros(dis,dis);
while err > 10^-8
    T_new= zeros(dis,dis);
    for j=2:dis-1
        for i=2:dis-1
            T_new(i,j)= lambda_y*T(i,j-1)+lambda_x*T(i-1,j)...
                +(1-2*lambda_x-2*lambda_y)*T(i,j)...
                +lambda_x*T(i+1,j)+lambda_y*T(i,j+1);
        end
    end
    % set the boundary conditions
    T_new(1,:)= x.*(b_x-x).^2;         % top side direchlet boundary condition
    T_new(dis,:)=(b_x-x).^2.*cos(pi*x/b_x);          % bottom side direchlet boundary condition
    gb_ax = (b_x-a_x)^2*cos(pi*a_x/b_x);
    fb_ax = a_x*(b_x-a_x)^2;
    T_new(:,1)= fliplr(gb_ax +(y-a_y)./(b_y-a_y)*(fb_ax-gb_ax));       % left side
    for i = 2:dis-1
        T_new(i,dis)=lambda_y*T(i,j-1)+2*lambda_x*T(i-1,j)...
            +(1-2*lambda_x-2*lambda_y)*T(i,j)+lambda_y*T(i,j+1);
    end
    
    T = T_new;
    time = time+dt;
    T_new = flipud(T_new);
    
    % error check:
    a = mean( mean(T_new));
    b = mean( mean(T_previous));
    err = abs( (a-b) / a * 100 );
    
    if ( mod(iter,1000) == 0)
        disp('Current iteration and error:')
        disp(iter)
        disp(err)
    end
    iter = 1+iter;
    T_previous = T_new;
    
end