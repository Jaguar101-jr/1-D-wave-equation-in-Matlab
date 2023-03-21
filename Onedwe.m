%% This code is written by HUSSEIN ABDUELHALEIM HUSSEIN MUHAMMED, BScH from UofK-Sudan,
%% Now at UPC (China), Aug. 2021 as part of his master thesis's project.

%% A Finite-Difference program to solve the 1-D wave equation in Cartesian domain

%% Please cite this code as: Hussein Abduelhaleim Hussein Muhammed (2022). Least-Squares
%% Reverse Time Migration in Pseudodepth Domain and Its Application.China University of Petroleum (East China), Master Thesis.
%% School of Geosciences, Dept. of Geophysics, Library press.
%% Thanks is to Dr. Haroon Stephen.

clear;

%% Wave Equation: Wtt = (c^2 * Wxx) + f
%% Prepare the movie file
    vidObj = VideoWriter('WE-1d.avi');
    open(vidObj);
%% Domain
% Space
Lx=10;
dx=0.1;
nx=fix(Lx/dx);
x=linspace(0, Lx, nx);
 
%Time
T=7;
 
%% Field variable
% Variables
wn=zeros(nx,1);
wnm1=wn; % w at time n-1
wnp1=wn; % w at time n+1
 
% Parameters
CFL=0.1; % CFL = c.dt/dx
c=1;
dt=CFL*dx/c;
 
%% activate/deactivate the Initial Conditions
wn(49:51)=[0.1 0.2 0.1];
wnp1(:)=wn(:);

%Time stepping Loop
t=0;
 
while(t < T)
    
    %%choose a boundary condition to activate
    % Reflecting Boundary Conditions
    %wn([1 end])=0;
    
    % Mur's Absorbing Boundary
    % wnp1(1)=wn(2) +((CFL-1)/(CFL+1))*(wnp1(2)-wn(1));
    % wnp1(end)=wn(end-1) + ((CFL-1)/(CFL+1)*(wnp1(end-1)-wn(end)));
    
    % Exact Solution
    t=t+dt;
    wnm1=wn; wn=wnp1; % Save current and previous arrays
    
    % activate/deactivate the source wavelet 'f'
    wn(50)=dt^2*20*sin(20*pi*t/T);                          % Ricker
    
    for i=2:nx-1
        wnp1(i) = 2*wn(i)- wnm1(i) ...
                   + CFL^2 * (wn(i+1) - 2*wn(i) + wn(i-1));
    end
    
    % Visualize the solution at selected steps
    clf;
    plot(x, wn);
    title(sprintf('time elapsed = %.2f' , t));
    axis([0 Lx -0.05 0.05]);
    shg; pause(0.01);
    xlabel('distance');
    ylabel('Amplitude');
    
 % Write each frame to the file
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
end


%% Close the file
close(vidObj);
