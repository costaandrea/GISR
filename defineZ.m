function [Z, Y] = defineZ(zi,xi,m,n,zs,a,dx)


for ii=1:n %for every x
    
    zslope = -a*xi(ii)+zi(zs);%equation of the slope
    
    zz(ii)=zslope;
   if zslope>=0 %if the slope is still in the domain (z>0)
     for jj=1:m %for every z
    
       if zi(jj)>=zslope %above the slope
          
         Z(jj,ii) = zi(jj)-zslope; %height above the slope
       
        elseif zi(jj)<zslope %below the slope (ground)
%            
          Z(jj,ii)=0;%NaN; 
           
       end
     end%for on z
     
   elseif zslope<0 %the slope is no longer in the domain
       
       for jj=1:m %for every z
           
          Z(jj,ii) = zi(jj)-zslope; %height above the slope
          
       end
   end
end%for on x 

for jj=1:n
    
    yi(jj) = (jj-1)*dx;
    
end
Y = repmat(yi,size(Z,1),1);
