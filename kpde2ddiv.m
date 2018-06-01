function [B1,B2]=kpde2ddiv(p,t)
%KPDE2DDIV Assembly of 2D divergence matrices with P1 finite element
%--------------------------------------------------------------------
%       [B1,B2]=kpde2ddiv(p,t) 
%
%  Input:
%         p : Node coordinates, np*2
%         t : Triangle vertices, nt*3  
%
%  Output:
%     B1,B2 : Divergence matrices, sparse np*np
%--------------------------------------------------------------------
% (c) J. Koko, LIMOS 2006-2016, koko@isima.fr
%--------------------------------------------------------------------
np=size(p,1);
% gradients of basis functions & triangle areas
[vol,g] = kpde2dgphiold(p,t);  
c=vol/3; 
% cell-array of gradient
D1={g{1}(:,1) g{2}(:,1) g{3}(:,1)};
D2={g{1}(:,2) g{2}(:,2) g{3}(:,2)};
% assembly of divergence matrices
B1=sparse(np,np);
B2=sparse(np,np);
for i=1:3
for j=1:3
    B1=B1+sparse(t(:,i),t(:,j),c.*D1{j},np,np);
    B2=B2+sparse(t(:,i),t(:,j),c.*D2{j},np,np);
end
end
