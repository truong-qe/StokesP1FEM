function [A,B,C]=kstok3dp1bmat(p,t,nu,alpha,vol,g)
%KSTOK3DP1BMAT 3D Stokes problem with  P1-Bubble/P1 element
% assembly of the matrix
%--------------------------------------------------------------------
%          A=kstok3dp1bmat(p,t,nu,alpha)
%          A=kstok3dp1bmat(p,t,nu,alpha,vol,g)
%          [A,B,C]=kstok3dp1bmat(p,t,nu,alpha),      
%          [A,B,C]=kstok3dp1bmat(p,t,nu,alpha,ar,g) 
%
%  Input:
%         p : Node coordinates, np*3
%         t : Tetrahedron vertices, nt*4  
%        nu : Viscosity coefficient, scalar or column vector nt*1
%     alpha : Mass coefficient, scalar or column vector nt*1
%       vol : Tetrahedrons volume
%         g : ={g1 g2 g3 g4} Cell-array of gradient of basis functions
%
%  Output:
%         ----------- 1 output argument -----------------------------
%         A : Stokes matrix, sparse (4*np)*(4*np)
%             A=[ K   0   0  -B1'
%                 O   K   0  -B2'
%                 0   0   K  -B3'
%                -B1 -B2 -B3 -C ]
%             K the stiffness matrix, [B1 B2 B3] the divergence matrix, 
%             C  the pressure matrix
%         ------------ 3 output arguments ----------------------------
%         A : Stiffness matrix K, np*np
%         B : B=[B1 B2] divergence matrix, np*(3*np)
%         C : Pressure stabilisation matrix, np*np
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2016, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1); Z=sparse(np,np);
% Gradient of basis functions
if (nargin == 4) [vol,g]=kpde3dgphi(p,t); end
% Bubble coefficients ( z and omega )
zt=(8/105)*alpha.*vol;
omega=(8192/945)*nu.*vol.*(sum(g{1}.^2,2)+sum(g{2}.^2,2)+sum(g{3}.^2,2)...
     +sum(g{1}.*g{2},2)+sum(g{1}.*g{3},2)+sum(g{2}.*g{3},2))+(8192/51975)*alpha.*vol;
c=(32/105)*vol;
% Stiffness, mass and divergence matrices
Ah=sparse(np,np); Bh=cell(1,3); [Bh{:}]=deal(Z); Ch=sparse(np,np);
for i=1:4
   for j=1:4
       Ah=Ah+sparse(t(:,i),t(:,j),nu*vol.*sum(g{i}.*g{j},2),np,np)...
            +sparse(t(:,i),t(:,j),alpha*vol/20,np,np)...
            -sparse(t(:,i),t(:,j),zt.*zt./omega,np,np);
       Ch=Ch-sparse(t(:,i),t(:,j),c.*c.*sum(g{i}.*g{j},2)./omega,np,np);
       for k=1:3
           Bh{k}=Bh{k}-sparse(t(:,i),t(:,j),vol.*g{j}(:,k)/4,np,np)...
                      -sparse(t(:,i),t(:,j),c.*g{i}(:,k).*zt./omega,np,np);  
       end
   end
   Ah=Ah+sparse(t(:,i),t(:,i),alpha*vol/20,np,np); 
end
% Final matrix
if (nargout == 1)
    A=[Ah Z  Z  Bh{1}'; Z  Ah Z  Bh{2}'; Z  Z  Ah Bh{3}'; Bh{1} Bh{2} Bh{3} Ch];
elseif (nargout == 3)
    A=Ah; B=[-Bh{1} -Bh{2} -Bh{3}]; C=-Ch;
end
 

