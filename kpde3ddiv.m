function [B1,B2,B3]=kpde3ddiv(p,t,vol,g)
%KPDE3DDIV Assembly of 3D divergence matrices with P1 finite element
%--------------------------------------------------------------------
%  [B1,B2,B3]=kpde3ddiv(p,t)
%  [B1,B2,B3]=kpde3ddiv(p,t,vol,g)
%
%  Input:
%         p : Node coordinates, np*3
%         t : Tetrahedron vertices, nt*4  
%       vol : elements volume, nt*1
%         g : {g1 g2 g3 g4} cell-array of gradient of basis functions, nt*3
%
%  Output:
%  B1,B2,B3 : Divergence matrices, sparse np*np
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2016, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1);  

% compute gradients of basis functions if necessary
if (nargin == 2)
   [vol,g] = kpde3dgphi(p,t);
elseif (nargin  ~= 7), error('Not enough or too many input arguments'), end
c=vol/4; 
% cell-array of gradient
D1={g{1}(:,1) g{2}(:,1) g{3}(:,1) g{4}(:,1)};
D2={g{1}(:,2) g{2}(:,2) g{3}(:,2) g{4}(:,2)};
D3={g{1}(:,3) g{2}(:,3) g{3}(:,3) g{4}(:,3)};
% assembly of divergence matrices
B1=sparse(np,np); B2=sparse(np,np); B3=sparse(np,np);
for i=1:4
for j=1:4
    B1=B1+sparse(t(:,i),t(:,j),c.*D1{j},np,np);
    B2=B2+sparse(t(:,i),t(:,j),c.*D2{j},np,np);
    B3=B3+sparse(t(:,i),t(:,j),c.*D3{j},np,np);
end 
end
