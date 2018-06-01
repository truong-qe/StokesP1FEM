function kstokshow(p,t,u1,u2,u3,pr,bcc)

uu=sqrt(u1.^2+u2.^2+u3.^2);

xb=[0 1 1 0 0]; xh=xb;
yb=[0 0 1 1 0]; yh=yb;
zb=[0 0 0 0 0]; zh=zb-1;
xg=[0 0 0 0 0]; xd=xg+1; 
yg=[0 1 1 0 0]; yd=yg;
zg=[0 0 -1 -1 0]; zd=zg;
xf=[0 1 1 0 0]; xb=xf;
yf=[0 0 0 0 0]; yb=yf+1;
zf=[0 0 -1 -1 0]; zb=zf;



close all
kpde3dshow(p,t,pr) 
hold on
plot3(xb,yb,zb,'k-',xh,yh,zh,'k-',xg,yg,zg,'k-',xd,yd,zd,'k-',xf,yf,zf,'k-',xb,yb,zb,'k-')
hold on
quiver3(p(bcc,1),p(bcc,2),p(bcc,3),u1(bcc),u2(bcc),u3(bcc),2.5)
xlabel('x'),ylabel('y'),zlabel('z'), title('Pressure and Velocity field')
colorbar east
axis tight, colormap hsv
view(3)
