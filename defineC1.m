function C1 = defineC1(zi,m,zc,n,m0,m1)

for i=1:n
    for j=1:m
        
        if zi(j)>zi(zc)
             C1(j,i)=m1;
         else
             C1(j,i)=m0;
        end
    end
end