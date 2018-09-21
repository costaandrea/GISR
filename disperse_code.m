clear all
close all
clc

%%%path to directory
cd('C:\Users\andre\Desktop\GISR_sims')%C:\Users\Kurt\Desktop\ALL')
%%%load data for temporal tuning of K0
load('U3.mat','M6')%C:\Users\Kurt\Desktop\ALL\U3.mat','M6');

fprintf(['=== START PROGRAM ===\r\r'])
fprintf(['Creating domain...\r\r'])
%tell me
dx = 100;
dz = 2;

days = 4*30 -12;

a=atan(0.03);

%free parameters
Cd = 0.01;
beta = 1;
alpha = 2*pi;

coeffiso = 50;
Kint=1e-5;

%%%define domain
xi = 0:dx:100e3;%250e3;
zi = 0:dz:1000;
m=length(zi);
n=length(xi);

zc=findnearest(zi,750);

%create grid
[x,z] = meshgrid(xi,zi);

%z where slope starts
zs = find(zi==max(zi)); 
%set slope
 for i = 1:findnearest(xi,zi(zs)*50) %if u want bottom  :n
  xind(i) = i;
  zind(i) = findnearest(zi, fix(-xi(i)/50 +zi(zs)));
%  
 end
 
 
%%%define initial condition (concentration)
sigx=50; %streak 100m
sigz=5; %2m

xc = min( findnearest( xi, fix( -zi(zc)*50 + zi(zs)*50)  +7500 ) );

fprintf(['Creating initial condition...\r\r'])

load('mean_profile.mat','Cbar','hgrd')
Cbar=Cbar(1:5:end); 
hgrd=hgrd(1:5:end)+zi(zc);

cz=zeros(1,m);
for k=1:numel(hgrd)
    cz(zi==hgrd(k))=Cbar(k);
end
cz = repmat(cz',1,n); 

normx = normpdf(xi,xi(xc),100);
cx = repmat(normx,m,1);

c0 = cx.*cz; %gaussian initial condition

c0 = c0 ./ sum(sum(c0)) *16.8; %[kg]

% %%%verify c0 position
% figure
%     contour(xi/1e3,zi-zi(zc),c0)
%     %imagesc(xi/1e3,zi-zi(zc),c0)
%     colorbar
%     xlabel('x [km]')
%     ylabel('z [m]')
%     title(num2str(0))
%     set(gca,'fontsize',14)
%     line(xi(xind)/1e3,zi(zind)-zi(zc),'color','k','linewidth',1)
%     line(xi(xind)/1e3,zi(zind)-zi(zc)+100,'color','r','linestyle','--','linewidth',1)

fprintf(['Defining physical variables...\r\r'])

%%%define buoyancy
N2=zeros(m,n);

m0=(5e-6 - 1e-6)/400;
m1=(1e-5 - 5e-6)/375;
cc1=5e-6;
cc0=1e-6;
 for i=1:m  
     if zi(i)>=zi(zc)
         N2(i,:) = (1e-5 - 5e-6)/375*(zi(i)-zi(zc)) + N2(zc-1,1);%5e-6;
     else
         N2(i,:) = (5e-6 - 1e-6)/400*zi(i) + 1e-6;
     end
 end
%figure;plot(log10(N2(:,end)),zi)
%set(gca,'YDir','Reverse')
N=sqrt(N2);
%figure;imagesc(N);colorbar


 for ii=1:m
        xind_inv(ii,:) = max(findnearest( xi, fix( -zi(ii)*50 + zi(zs)*50) ));  
 end
    
%%%isopycanl diffusivity
K1=coeffiso*ones(m,n);
%normf=0.01;
%K1=getK1_old(xi,zi,zs,xind_inv,normf,m,n);

%%%make matrices of z, c0, c1
[Zm, Ym]=defineZ(zi,xi,m,n,zs,a,dx);  
Cs=defineC(zi,m,zc,n,cc0,cc1); 
Ms=defineC1(zi,m,zc,n,m0,m1); 

mask=ones(m,n);
mask(Zm==0)=0;

K1=K1.*mask;
N=N.*mask;

%%%advection time
T = 3600*24*days; %[s]

b = waitbar(0,'Initializing waitbar...');

cits(:,:,1)=c0;%initial condition
k=1;

del=0.5; Th=3600; %to have an output from advection scheme every hour

fprintf(['DISPERSING!\r\r'])

%figure
for timeind=1:fix(T/3600)

%%%flow velocity at this time    
U = nthroot(M6.U3(16*24 +12*24 +timeind),3)/100;
 
%%%diapycnal diffusivity
K3=getKz(alpha,beta,Cd,U,Kint,Zm,N).*mask;
%K3=getKz_bench(Cd,U,beta,Kint,N,Zm).*mask;
K3(isnan(K3))=0; %due to 1/n with N=0 below the slope
%figure;imagesc(K3);colorbar
%return

%%%diapycnal velocity
w=getW(alpha,beta,Cd,U,Kint,Zm,N,Ms).*mask;
%w=getW_bench(m,n,a,xi,zi,zs,Cd,U,beta,N,Zm,Cs,Ms,Kint).*mask;
w(Zm==0)=0; %no-flux across the slope
%figure;imagesc(w);colorbar


%%%isopycnal velocity (numerical integration)
u= diff(w,1)/dz*dx; %from contunuity equation
u= [u(1,:); u].*mask;
%u(Zm==0)=0; %no-flux across the slope %tested this. not needed.
%figure;imagesc(u);colorbar


%%%stability criteria
D1=max(max(K1));
D3=max(max(K3));
dt = 1/( 2*D1/dx^2 + 2*D3/dz^2 + max(max(abs(u)))/dx + max(max(abs(w)))/dz );

cycl = fix(T/3600);

%%%call advection scheme
cT  = diff_2D_ULT_vel(Th,timeind,cycl,dt,dx,dz,xi,zi,m,n,K1,K3,u,w,c0,xi,zi,xc,zc,sigz,zs,xind,zind,xind_inv);
c0 = cT;
%imagesc(flipud(c0))      
%colorbar



TotC(timeind) = sum(sum(c0));
    
if rem(timeind,24)==0  
      k=k+1;
     cits(:,:,k)=c0;
end
 
waitbar(1/fix(T/3600)*timeind,b,[num2str(timeind),' / ',num2str(cycl)])
end
close(b)
fprintf(['Done.\r\r'])

% figure
%     contour(xi/1e3,zi-zi(zc),c0)
%     colorbar
%     xlabel('x [km]')
%     ylabel('z [m]')
%     title(num2str(0))
%     set(gca,'fontsize',14)
%     line(xi(xind)/1e3,zi(zind)-zi(zc),'color','k','linewidth',1)


%%%save stuff
cd RUNS/
% name=['Days_',num2str(days),'_Cd',num2str(Cd),'_beta',num2str(beta),'_alp',num2str(alpha),...
%     '_dx',num2str(dx),'_dz',num2str(dz),'_Kint',Kint,'_Kiso',coeffiso,'.mat'];
name=['Cd',num2str(Cd),'B',num2str(beta),'A',...
    num2str(floor(alpha*100)/100),...
    'DX',num2str(dx),'DZ',num2str(dz),'KZ',...
    num2str(abs(floor(log10(Kint))))...
    ,'KY',num2str(coeffiso),'.mat'];

fprintf(['Savin` it...\r\r'])
save(name);

fprintf(['The file name is:\r',name(1:end-4),'\r\r'])

fprintf(['=== END PROGRAM ===\r'])
fprintf(['Thank you for running me!\r\r'])
