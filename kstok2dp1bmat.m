function [A,B,C]=kstok2dp1bmat(p,t,nu,alpha,ar,g)
%STOK2DP1BMAT Stokes problem with  P1-Bubble/P1 element
% assembly of the matrix
%--------------------------------------------------------------------
%          A=kstok2dp1bmat(p,t,nu,alpha),      
%          A=kstok2dp1bmat(p,t,nu,alpha,ar,g) 
%          [A,B,C]=kstok2dp1bmat(p,t,nu,alpha),      
%          [A,B,C]=kstok2dp1bmat(p,t,nu,alpha,ar,g) 
%
%  Input:
%         p : Node coordinates, np*2
%         t : Triangle vertices, nt*3  
%        nu : Viscosity coefficient, scalar or column vector nt*1
%     alpha : Mass coefficient, scalar or column vector nt*1
%        ar : Triangles area
%         g : ={g1 g2 g3} Cell-array of gradient of basis functions
%
%  Output:
%         ----------- 1 output argument -----------------------------
%         A : 1-output argument, Stokes matrix, sparse (3*np)*(3*np)
%             A=[ K   0  -B1'
%                 O   K  -B2'
%                -B1 -B2 -C ]
%             K the stiffness matrix, [B1 B2] the divergence matrix, 
%             C  the pressure matrix
%         ------------ 3 output arguments ----------------------------
%         A : Stiffness matrix K, np*np
%         B : B=[B1 B2] divergence matrix, np*(2*np)
%         C : Pressure stabilisation matrix, np*np
%         
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2016, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); Z=sparse(np,np);
% Triangles area and gradient of basis functions
if (nargin == 4) [ar,g]=kpde2dgphiold(p,t);  end   
% Bubble coefficients z and omega
zt=(3/20)*alpha.*ar;
omega=(81/10)*nu.*ar.*(sum(g{1}.^2,2)+sum(g{2}.^2,2)+sum(g{1}.*g{2},2))...
     +(81/280)*alpha.*ar;
c=(9/10)*ar;
% Mass & stiffness matrices
Ah=Z;  Bh=cell(1,2); [Bh{:}]=deal(Z); Ch=sparse(np,np);
for i=1:3
    for j=1:3       
       Ah=Ah+sparse(t(:,i),t(:,j),nu.*ar.*sum(g{i}.*g{j},2),np,np)...
            +sparse(t(:,i),t(:,j),alpha.*ar/12,np,np)...
            -sparse(t(:,i),t(:,j),zt.*zt./omega,np,np);
       Ch=Ch-sparse(t(:,i),t(:,j),c.*c.*sum(g{i}.*g{j},2)./omega,np,np);
       for k=1:2
           Bh{k}=Bh{k}-sparse(t(:,i),t(:,j),ar.*g{j}(:,k)/3,np,np)...
                      -sparse(t(:,i),t(:,j),c.*g{i}(:,k).*zt./omega,np,np);
       end
    end
    Ah=Ah+sparse(t(:,i),t(:,i),alpha.*ar/12,np,np);
end
% 
% Output
%
if (nargout == 1)
    A=[Ah   Z   Bh{1}'; Z    Ah  Bh{2}'; Bh{1}   Bh{2}  Ch];
elseif (nargout > 1)
    A=Ah; B=[-Bh{1} -Bh{2}]; C=-Ch;
end
    

