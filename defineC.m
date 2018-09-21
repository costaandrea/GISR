function C = defineC(zi,m,zc,n,c0,c1)

for i=1:n
    for j=1:m
        
        if zi(j)>zi(zc)
             C(j,i)=c1;
         else
             C(j,i)=c0;
        end
    end
end