function y = diff_2D_ULT_vel(T,timeind,cycl,dt,dx,dz,x,z,m,n,K1,K3,u,v,c0,xi,zi,xc,zc,sigx,zs,xind,zind,xind_inv)
% y=diff_2D(x,z,T,K1,K3,u,c0);
% Created on 8/4/02; Ledwell; for SFTRE scenarios
% Advection Diffusion in an infinite medium - 2-D
% with diffusivities K1 and K3 depending on x and z
% x and z are evenly spaced grids
% and a zonal velocity u depending on z
% for time T on even grid x,z
% initial condition vector c0
% gradient = 0 at the boundaries

% Solves the diffusion equation
% using an explicit algorithm
% Solution y(x,z) is returned at the true final time

% Stability requires dx^2>2*K1max*dt where K1max = max(K1)
% and dz^2>2*K3max*dt where K3max=max(K3);


%%%NB Adding a fake line/row
kz=[K3(2,:); K3; K3(m-2,:)]; %matrix of K3
kx=[K1(:,2) K1 K1(:,n-2)]; %matrix of K1
p=[u(:,2) u u(:,n-2)]; %matrix of u

q=[v(2,:); v; v(m-2,:)];%matrix of v


r1=dt/(2*dx); %
r3=dt/(2*dz); %

y=c0; %solution equal to initial condition

%%%NB Adding a fake line/row
Z=zeros(m,1); 
Z3=zeros(1,n); 

%%%operators
% wu1 = r1*[Z (p(:,1:n-2)+p(:,2:n-1))/2 + (kx(:,1:n-2)+kx(:,2:n-1))/dx Z];
% wu3 = r3*[Z3; (kz(1:m-2,:) + kz(2:m-1,:))/dz; Z3];
% 
% wc1 = r1*[Z (p(:,1:n-2)-p(:,3:n))/2 - (kx(:,1:n-2)+2*kx(:,2:n-1)+kx(:,3:n))/dx Z];
% wc3 = r3*[Z3; -(kz(1:m-2,:)+2*kz(2:m-1,:)+kz(3:m,:))/dz; Z3];
% 
% wd1 = r1*[Z -(p(:,2:n-1)+p(:,3:n))/2 + (kx(:,2:n-1)+kx(:,3:n))/dx Z];
% wd3 = r3*[Z3; (kz(2:m-1,:)+kz(3:m,:))/dz; Z3];

 
wu1 = r1*[Z (p(:,1:n-2)+p(:,2:n-1))/2 + (kx(:,1:n-2)+kx(:,2:n-1))/dx Z];
wu3 = r3*[Z3;  (q(1:m-2,:)+q(2:m-1,:))/2 + (kz(1:m-2,:) + kz(2:m-1,:))/dz; Z3];

wc1 = r1*[Z (p(:,1:n-2)-p(:,3:n))/2 - (kx(:,1:n-2)+2*kx(:,2:n-1)+kx(:,3:n))/dx Z];
wc3 = r3*[Z3; (q(1:m-2,:)-q(3:m,:))/2 - (kz(1:m-2,:)+2*kz(2:m-1,:)+kz(3:m,:))/dz; Z3];

wd1 = r1*[Z -(p(:,2:n-1)+p(:,3:n))/2 + (kx(:,2:n-1)+kx(:,3:n))/dx Z];
wd3 = r3*[Z3; -(q(2:m-1,:)+q(3:m,:))/2 + (kz(2:m-1,:)+kz(3:m,:))/dz; Z3];

%integration

 %slope_ind = [xind,zind];
 
% slope = ones(m,n);
% for i=1:n
%     slope(i,zind:n)=0;
% end

%%%time dependence of diffusivities from U3
% load('/Users/andrea/Desktop/ANDREA@WHOI/U3.mat');
% U3 = squeeze(M6.U3);

N=1.7e-3;

Nt = fix(3600/dt);

k=0;
kk=0;

for j=1:Nt-1
   
 %%%every hour change K3 calling K3maker.m   
%  if rem(fix(j*dt/0.041),0.041)==0
%      Uind = fix(j*dt/0.041);
%      K3 = K3maker(Uind,N,zi,zs,U3);
%  end
%    
%     if dt*j<timeind/24
%        K3 = K3maker(timeind,N,zs,U3,m,n,xi,zi,xind_inv);
%        
%     else 
%         timeind = timeind+1;
%         K3 = K3maker(timeind,N,zs,U3,m,n,xi,zi,xind_inv);
%     end

 %plot(zi,y(:,10),'-bo')
 y = y + wu1.*[Z y(:,1:n-1)] + (wc1+wc3).*y + wd1.*[y(:,2:n) Z] + ...
     wu3.*[Z3; y(1:m-1,:)] + wd3.*[y(2:m,:); Z3];
 
     %%% no flux b.c. %JIM code
     %y(:,1)=y(:,3);  %%%%why not y(:,2)? etc.%%%%because of the fake line%%%
     %y(1,:)=y(3,:); %vertical left wall
     %y(:,n)=y(:,n-2); %rigth margin
     %y(m,:)=y(m-2,:); %top
     
     %%%put b.c. on the slope
     for ii=1:numel(zind)
         y(zind(ii),ii) = y(min(zind(ii)+2,m),ii); 
     end  

%     if rem(j,1000)==0 %u don wanna blow up the memory
%        k=k+1;
%         
%         g(k,:) = gradient(y(zc,:));
%         yp(k,:,:) = y;
%         if sum(sum(y.*VC))>0.001
%             kk=kk+1;
%             yp_ns(kk,:,:)=y; %near slope
%         end
        
%     end

end
