clc;
clear;
tic;
% Crank-Nicholson method is unconditionally stable
nx=41; % Number of grid points in X-direction
ny=41; % Number of grid points in Y-direction
dt=0.0001; % Time-step 
T=0.16;
dx=1/(nx-1); % Width of one grid point
dy=1/(ny-1); % Width of one grid point
t=0:dt:T;
[row_t,column_t]=size(t);% Storing the size of time step to create the time evolution plot
w3=zeros(column_t,1); % Column vector to store the temperature value of x=y=0.4 at each time step from 0 to 0.16
% Grid points spacing in X direction
x=0:dx:1;
% Grid points spacing in Y direction
y=0:dy:1;
% a is the constant in the stencil = (2*h*h)/delta t
a=(2*dx*dx)/dt;
% Matrix to store the n+1 time step values of the Temperature
w1= zeros((nx-2)*(ny-2),1);
% Matrix to store the new step value and use it for the next iteration
w5= zeros((nx-2)*(ny-2),1); 
% Matrix to store the boundary conditions
w= zeros(nx,ny); 
w(1,:)=1;
w(nx,:)=0;
w(:,1)=1-y;
w(:,nx)=1-y.*y;

% RHS is the right side constants in stencil 
% This will differ from the boundary matrix as the RHS matrix contains both
% values of the boundary at n+1 and n time step
rhs=zeros((nx-2),(ny-2));


% Initialization of RHS values begin
% To store the left side boundaries
for j=1
    for i=1:(nx-2)
        if i==j
            rhs(i,j)=(2*w(i,j+1))+(2*w(i+1,j));
        else
            rhs(i,j)=(2*w(i+1,j));
        end
    end
end

%To store the top boundary values
for i=1
    for j=1:(nx-2)
        if i~= j
            rhs(i,j)=(2*w(i,j+1));
        end
    end
end

%To store the bottom boundary values
rhs(end,2:end-1)=0;

%To store the right side boundary values
for j=(nx-2)
    for i=1:(nx-2)
        if i==1
            rhs(i,j)=(2*w(1,1))+(2*w(i+1,j+2));
        else
            rhs(i,j)=(2*w(i+1,j+2));
        end
    end
end

% Initialization of boundary value ends
   
% To create the LHS Tri-diagonal matrix 
T1=diag((a+4)*ones(1,(nx-2)*(nx-2)),0)+diag(-1*ones(1,(nx-2)*(nx-2)-1),1)+diag(-1*ones(1,(nx-2)*(nx-2)-1),-1);
n=nx-2; % New variables initialized for the purpose of if loops inside for loops
n1=nx-2; % New variables initialized for the purpose of if loops inside for loops
for i=1:(nx-2)*(nx-2) % for loop to modify the required elements in the tri-diagonal matrix
    for j=1:(nx-2)*(nx-2)
        if n1 ~= (nx-2)*(nx-2)  % This condition is checked, so that when n1 reaches nx*nx, it will 
                        % not execute the below if condition 
            if j==n1 && i==n1 % if condition to assignment the certain upper-diagonal and 
                              % lower diagonal elements to zero according
                              % to the matrix i created for 3*3 grid points 
            T1(i+1,j)=0;
            T1(i,j+1)=0;
            n1=n1+(nx-2);
            end
        end
        if j==n+1  % if condition to assign the other minor diagonals of the matrix elements
        T1(i,j)=-1; % These elements which are assigned here are symmetric based on the 
        T1(j,i)=-1; % example matrix i created during descritization of the heat equation
        end
    end
    n=n+1;
end

% To create the RHS tri-diagonal matrix
T2=diag((a-4)*ones(1,(nx-2)*(nx-2)),0)+diag(1*ones(1,(nx-2)*(nx-2)-1),1)+diag(1*ones(1,(nx-2)*(nx-2)-1),-1);
n=(nx-2); %N is again initialized to nx as it was modified in the previous tri-diagonal creation
      % New variables initialized for the purpose of if loops inside for loops
n2=(nx-2); % New variables initialized for the purpose of if loops inside for loops
for i=1:(nx-2)*(nx-2)
    for j=1:(nx-2)*(nx-2)
        if n2 ~= (nx-2)*(nx-2)  % This condition is checked, so that when n2 reaches nx*nx, it will 
                        % not execute the below if condition 
            if j==n2 && i==n2
            T2(i+1,j)=0;
            T2(i,j+1)=0;
            n2=n2+(nx-2);
            end
        end
        if j==n+1 % if condition to assign the other minor diagonals of the matrix elements
            T2(i,j)=1; % These elements which are assigned here are symmetric based on the 
            T2(j,i)=1; % example matrix I created during descritization of the heat equation
        end
    end
    n=n+1;
end

rhs=reshape(rhs,(nx-2)*(nx-2),1); % To reshape the rhs matrix into a column vector

%Solving LSE
for t1=0:dt:T
     
    
    w1=T1\(T2*w5+rhs);
    
    
    if t1==0.01
        w2=reshape(w1,(nx-2),(ny-2));
        w(2:end-1,2:end-1)=w2;
        figure('Name','Temp Contour Plot at t = 0.01 by Crank-Nicolson','NumberTitle','off')
        contourf(x,y,w);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.01','FontSize',24,'FontWeight','bold')
    end
     if t1==0.02
        w2=reshape(w1,(nx-2),(ny-2));
        w(2:end-1,2:end-1)=w2;
        figure('Name','Temp Contour Plot at t = 0.02 by Crank-Nicolson','NumberTitle','off')
        contourf(x,y,w);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.02','FontSize',24,'FontWeight','bold')
     end
     if t1==0.04
        w2=reshape(w1,(nx-2),(ny-2));
        w(2:end-1,2:end-1)=w2;
        figure('Name','Temp Contour Plot at t = 0.04 by Crank-Nicolson','NumberTitle','off')
        contourf(x,y,w);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.04','FontSize',24,'FontWeight','bold')
     end
     if t1==0.08
        w2=reshape(w1,(nx-2),(ny-2));
        w(2:end-1,2:end-1)=w2;
        figure('Name','Temp Contour Plot at t = 0.08 by Crank-Nicolson','NumberTitle','off')
        contourf(x,y,w);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.08','FontSize',24,'FontWeight','bold')
     end
     w4=reshape(w1,(nx-2),(ny-2));
     w(2:end-1,2:end-1)=w4;
     i=fix((t1*10000)+1); % To eliminate decimals in i, if at all anything is there
     w3(i,1)=w(17,17); % corresponding to x=y=0.4 ,w is (1tth row and 17th column )
     if t1==0.16
         figure('Name',' Time evolution of Temperature at x=y=0.4 by Crank-Nicolson with stable delta t=0.0001','NumberTitle','off')
         plot(t,w3,'Linewidth',2)
         xlabel(' Time ','FontSize',18,'FontWeight','bold')
         ylabel('Temperature','FontSize',18,'FontWeight','bold')
         axis([0 0.16 0 1])
         title(' Time evolution of the temperature @ x=y=0.4 with delta t=0.0001','FontSize',24,'FontWeight','bold')
     end
     w5=w1;
end

w1=reshape(w1,(nx-2),(ny-2)); % To reshape the column vector into matrix
w(2:end-1,2:end-1)=w1; % To copy the calculated values into the boundary matrix

figure('Name','Temp Contour Plot at t = 0.16 by Crank-Nicolson','NumberTitle','off')
contourf(x,y,w);colorbar
xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.16','FontSize',24,'FontWeight','bold')

figure('Name','Temp vertical profile at x = 0.4and t = 0.16 by Crank-Nicolson','NumberTitle','off')
plot(y,w(:,17),'Linewidth',2);
xlabel('Y axis','FontSize',18,'FontWeight','bold')
ylabel('Temperature value','FontSize',18,'FontWeight','bold')
title('Vertical Temperature profile at x=0.4 , t = 0.16 , dt=0.0001','FontSize',24,'FontWeight','bold')
toc



