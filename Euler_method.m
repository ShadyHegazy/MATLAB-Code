
clc;
clear;
%Explicit Euler method
nx=41; % Number of grid points in X-direction
ny=41; % Number of grid points in Y-direction
dt=0.0001; % Time-step 
dx=1/(nx-1); % Width of one grid point
dy=1/(ny-1); % Width of one grid point
% To check stability of the scheme. if not stable, the problem quits
if dt<dx*dx/4
    fprintf('Present time step %f is stable \n',dt);
else
    error('Present time step %f is not stable \n Please give a new time step \n',dt);
end
t=0:dt:0.16;
[row_t,column_t]=size(t);
x=0:dx:1;
y=0:dy:1;
w= zeros(nx,ny);  % Matrix to store the n+1 time step values of the Temperature
wn= zeros(nx,ny); % Matrix to store the n time step values of the temperature
wn1=zeros(column_t,1);
a=dt/(dx*dy); % Constant A = delta T/ delta h^2
% Matrix to store the iteration results
wn(1,:)=1;
wn(ny,:)=0;
wn(:,1)=1-y;
wn(:,nx)=1-y.*y;
w(1,:)=1;
w(nx,:)=0;
w(:,1)=1-y;
w(:,nx)=1-y.*y;
tic;
for t1=0:dt:0.16 % For loop to calculate the temperature for every time step till the Maximum range
   
    for i=2:nx-1
        for  j=2:ny-1
         % For exach node points the n+1 time step temperature is
         % calculated from n time steps and stored in the matrix
         w(i,j)=((1-4*a)*wn(i,j))+a*(wn(i+1,j)+wn(i-1,j)+wn(i,j+1)+wn(i,j-1));
         % To generate the plot for Temperature evolution with respect to time
         % at x,y=0.4 with stable delta T
         if i==16 && j==16
            wn1(fix((t1*10000)+1),1)=w(i,j);
         end
        end
    end  
    wn=w; % Before the next iteration for next time step, copy the n+1 
          % time step temperature to n timestep temperature matrix and use
          % it for the next iteration
    if t1==0.01
        figure('Name',' Contour graphs for Euler Method at t=0.01','NumberTitle','off')
        contourf(x,y,wn);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.01','FontSize',24,'FontWeight','bold')
    end
    if t1==0.02
        figure('Name',' Contour graphs for Euler Method at t=0.02','NumberTitle','off')
        contourf(x,y,wn);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.02','FontSize',24,'FontWeight','bold')
    end
    if t1==0.04
        figure('Name',' Contour graphs for Euler Method at t=0.04','NumberTitle','off')
        contourf(x,y,wn);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.04','FontSize',24,'FontWeight','bold')
    end
    if t1==0.08
        figure('Name',' Contour graphs for Euler Method at t=0.08','NumberTitle','off')
        contourf(x,y,wn);colorbar
        xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
        ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
        title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.08','FontSize',24,'FontWeight','bold')
    end
end
  % contour graph for the Time step of 0.0001
  figure('Name',' Contour graphs for Euler Method at t=0.16','NumberTitle','off')
  contourf(x,y,wn);colorbar
  xlabel(' Spatial co-ordinate (x)','FontSize',18,'FontWeight','bold')
  ylabel('Spatial co-ordinate (y)','FontSize',18,'FontWeight','bold')
  title('2-D Diffusion Heat equation @ dt=0.0001 & t=0.16','FontSize',24,'FontWeight','bold')
  % Time evolution of Temperature at x=y=0.4 with stable delta T
  figure('Name',' Time evolution of Temperature at x=y=0.4 by Euler Method with stable delta t=0.0001','NumberTitle','off')
  plot(t,wn1,'Linewidth',2);
  xlabel( 'Time ','FontSize',18,'FontWeight','bold')
  ylabel(' Temperature ','FontSize',18,'FontWeight','bold')
  title(' Time evolution of Temperature at x=y=0.4 by Euler Method with stable delta t=0.0001','FontSize',24,'FontWeight','bold')
  % Vertical Temperature Profile at x=0.4 and t=0.16 with stable delta T
  figure('Name',' Vertical Temperature Profile at x=0.4 and t=0.16 with stable delta T','NumberTitle','off')
  plot(y,wn(:,17),'Linewidth',2);
  xlabel( 'Time ','FontSize',18,'FontWeight','bold')
  ylabel(' Temperature ','FontSize',18,'FontWeight','bold')
  title(' Vertical Temperature Profile at x=0.4 and t=0.16 with stable delta T','FontSize',24,'FontWeight','bold')
  toc

   
