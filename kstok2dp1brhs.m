function [b,bp]=kstok2dp1brhs(p,t,f1,f2,nu,alpha,ar,g)
%KSTOK2DP1BRHS  Stokes problem with  P1-Bubble/P1 element: 
% assembly of the right-hand side
%--------------------------------------------------------------------
%          b=kstok2dp1brhs(p,t,f1,f2,nu,alpha)
%          b=kstok2dp1brhs(p,t,f1,f2,nu,alpha,ar,g)
%          [b,bp]=kstok2dp1brhs(p,t,f1,f2,nu,alpha)
%          [b,bp]=kstok2dp1brhs(p,t,f1,f2,nu,alpha,ar,g)
%
%  Input:
%         p : Node coordinates, np*2
%         t : Triangle vertices, nt*3
%     f1,f2 : Source terms, scalar or column vector nt*1 or np*1
%        nu : Viscosity coefficient, scalar or column vector nt*1
%     alpha : Mass coefficient, scalar or column vector nt*1
%        ar : Triangles area
%         g : ={g1 g2 g3} Cell-array of gradient of basis functions
%
%  Output:
%         ------------ 1 output argument ----------------
%         b : Right-hand side, column vector (3*np)
%         ------------ 2 output argument ----------------
%         b : Velocity right-hand side, clolumn vector (2*np)
%        bp : Pressure right-hand side, column vector np
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2016, koko@isima.fr
%-------------------------------------------------------------------- 
np=size(p,1);  
% (f1,f2) at the center of triangles
if (length(f1)==np), f1t=sum(f1(t),2)/3;  else f1t=f1; end
if (length(f2)==np), f2t=sum(f2(t),2)/3;  else f2t=f2; end
% Triangles area and gradient of basis functions
if (nargin == 6) [ar,g]=kpde2dgphiold(p,t);  end
% z and omega  
zt=(3/20)*alpha.*ar;
omega=(81/10)*nu.*ar.*(sum(g{1}.^2,2)+sum(g{2}.^2,2)+sum(g{1}.*g{2},2))...
     +(81/280)*alpha.*ar;
c=(9/10)*ar; ft={f1t f2t}; 
% Assembly of the right-hand side
bb=sparse(np,1); bh=cell(2,1);[bh{:}]=deal(sparse(np,1));   % 
for i=1:3
for k=1:2
    bh{k}=bh{k}+sparse(t(:,i),1,(1/3)*ft{k}.*ar-c.*ft{k}.*zt./omega,np,1);
    bb=bb-sparse(t(:,i),1,c.*c.*g{i}(:,k).*ft{k}./omega,np,1);
end
end
% Right-hand side
if (nargout == 1) 
    b=[full(cell2mat(bh)); full(bb)];
elseif (nargout == 2)
    b=full(cell2mat(bh)); bp=full(bb);
end