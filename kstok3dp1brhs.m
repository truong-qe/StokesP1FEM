function [b,bp]=kstok3dp1brhs(p,t,f1,f2,f3,nu,alpha,vol,g)
%KSTOK3DP1BRHS  3D Stokes problem with  P1-Bubble/P1 element: 
% assembly of the right-hand side
%--------------------------------------------------------------------
%          b=kstok3dp1brhs(p,t,f1,f2,f3,nu,alpha)
%          b=kstok3dp1brhs(p,t,f1,f2,f3,nu,alpha,vol,g)
%          [b,bp]=kstok2dp1brhs(p,t,f1,f2,nu,alpha)
%          [b,bp]=kstok2dp1brhs(p,t,f1,f2,nu,alpha,ar,g)
%
%  Input:
%         p : Node coordinates, np*3
%         t : Tetrahedron vertices, nt*4
%  f1,f2,f3 : Source terms, scalar or column vector nt*1 or np*1
%        nu : Viscosity coefficient, scalar or column vector nt*1
%     alpha : Mass coefficient, scalar or column vector nt*1
%       vol : Tetrahedrons volume
%         g : ={g1 g2 g3 g4} Cell-array of gradient of basis functions
%
%  Output:
%         ------------ 1 output argument ----------------
%         b : Right-hand side, column vector (4*np)
%         ------------ 2 output argument ----------------
%         b : Velocity right-hand side, clolumn vector (3*np)
%        bp : Pressure right-hand side, column vector np
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2016, koko@isima.fr
%-------------------------------------------------------------------- 
np=size(p,1);
% (f1,f2) at the center of triangles
if (length(f1)==np), f1t=sum(f1(t),2)/4;  else f1t=f1; end
if (length(f2)==np), f2t=sum(f2(t),2)/4;  else f2t=f2; end
if (length(f3)==np), f3t=sum(f3(t),2)/4;  else f3t=f3; end
% Triangles area
if (nargin == 7) [vol,g]=kpde3dgphi(p,t); end
% z and omega  
zt=(8/105)*alpha.*vol; 
omega=(8192/945)*nu.*vol.*(sum(g{1}.^2,2)+sum(g{2}.^2,2)+sum(g{3}.^2,2)+sum(g{1}.*g{2},2)...
     +sum(g{1}.*g{3},2)+sum(g{2}.*g{3},2))+(8192/51975)*alpha.*vol;
c=(32/105)*vol; 
% Assembly of the right-hand side
ft={f1t f2t f3t};
bb=sparse(np,1);  bh=cell(3,1); [bh{:}]=deal(sparse(np,1));  
for i=1:4
for k=1:3
    bh{k}=bh{k}+sparse(t(:,i),1,(1/4)*ft{k}.*vol-c.*ft{k}.*zt./omega,np,1);
    bb=bb+sparse(t(:,i),1,c.*c.*g{i}(:,k).*ft{k}./omega,np,1);
end
end
% Output
if (nargout == 1)     b=[full(cell2mat(bh)); full(bb)];
elseif (nargout == 2) b=full(cell2mat(bh)); bp=full(bb); end