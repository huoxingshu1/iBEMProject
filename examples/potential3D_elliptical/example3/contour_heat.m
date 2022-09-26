clear
clc
x_num=;
y_num=;

x0=[0 0 0];
leng=0.6;


[A , B]=input1();

t=1;
for m=1:x_num
    for n=1:y_num
    xp(m,n) = B(t,2);      
    yp(m,n) = B(t,3);
    zp(m,n)= A(t);
    t=t+1;
    end
end

contourf(xp,yp,zp,30,'--k');
%contourf(xp,yp,zp,140,'edgecolor','none');

axis equal
shading flat
colorbar()
