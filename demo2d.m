%
% 2D Stokes Problem with P1-Bubble/P1 FEM
%
clear all, close all
%
alpha=0;
nu=0.01;
%
% choose method 1: LDL factorization, 2: Uzawa Conjugate Gradient
meth=2;   

% load mesh
load cyl2d1730 p t ibc1e ibc1i
np=size(p,1); nt=size(t,1);
fprintf('Mesh : nodes=%4d  triangles=%4d \n',np,nt)

%
% Assemble matrices and right-hand side
if (meth == 1)
    A=kstok2dp1bmat(p,t,nu,alpha);
    b=kstok2dp1brhs(p,t,zeros(np,1),zeros(np,1),nu,alpha);
elseif (meth == 2)
    [K,B,C]=kstok2dp1bmat(p,t,nu,alpha); 
    [b,bp]=kstok2dp1brhs(p,t,zeros(np,1),zeros(np,1),nu,alpha);
    A1=K; A2=A1;
    
    % Preconditioner
    K0=kpde2dstf(p,t,1); M0=kpde2dmss(p,t,1);
    DM0=spdiags(M0,0); DM0=spdiags(DM0,0,np,np);
    ar=kpde2dgphi(p,t);
    hx2=max(ar); Rnh=alpha*hx2/nu;
    Kh=K0+alpha*C; KM=nu*K0; 
    if (Rnh <= 1)
       KM=KM+alpha*DM0;
    else
       KM=KM+alpha*speye(np);
    end
end

%
% Incorporate the boundary conditions
PenBC=10^15;
if (meth == 1)
   ibc1=union(ibc1e,ibc1i); 
   ibc2=np+ibc1; ibcp=2*np+1;
   u1bc=0.3*4*p(ibc1e,2).*(0.41-p(ibc1e,2))/(0.41^2);
   A(ibc1,ibc1)=A(ibc1,ibc1)+PenBC*speye(length(ibc1));
   A(ibc2,ibc2)=A(ibc2,ibc2)+PenBC*speye(length(ibc2));
   A(ibcp,ibcp)=A(ibcp,ibcp)+PenBC;
   b(ibc1e)=PenBC*u1bc; b(ibc2)=0;
elseif (meth == 2)
   ibc1=union(ibc1e,ibc1i); 
   ibc2=ibc1; ibcp=1;
   u1bc=0.3*4*p(ibc1e,2).*(0.41-p(ibc1e,2))/(0.41^2);
   A1(ibc1,ibc1)=A1(ibc1,ibc1)+PenBC*speye(length(ibc1)); 
   A2(ibc2,ibc2)=A2(ibc2,ibc2)+PenBC*speye(length(ibc2)); 
   Kh(ibcp,ibcp)=Kh(ibcp,ibcp)+PenBC*speye(length(ibcp));
   b(ibc1e)=PenBC*u1bc;  bp(ibcp)=0;
end

tic
%
% Compute the solution (LDL' factorization)
if (meth == 1)
   [L,D,s]=ldl(A,'vector');
   warning('off')
   u=zeros(3*np,1);  u(s) = L'\(D\(L\(b(s))));
   u1=u(1:np); u2=u(np+1:2*np); pr=u(2*np+1:3*np);
   fprintf('LDL factorization CPU=%12.6f \n',toc)
elseif (meth == 2)
    [R1,p1,s1]=chol(A1,'lower','vector');
    [R2,p2,s2]=chol(A2,'lower','vector');
    [Rh,pp,sh]=chol(Kh,'lower','vector');
    
    % Data structure
    R={R1,R2}; s={s1 s2}; ibcd={ibc1 ibc2};  
    b={b(1:np) b(np+1:2*np)};
    [u,pr,iter]=kstokcg(R,Rh,B,C,KM,s,sh,b,bp,ibcd,ibcp,10^-6);
    u1=u(1:np); u2=u(np+1:2*np);
    fprintf('Uzawa Conjugate Gradient : iter=%5d CPU=%12.6f\n',iter,toc)
 end

%
% Compute the streamlines
K=kpde2dstf(p,t,1);
[B1,B2]=kpde2ddiv(p,t);
b=B2'*u1-B1'*u2;
K(1,1)=PenBC;
w=K\b;

%
% Outputs
figure('Name','Mesh')
   triplot(t,p(:,1),p(:,2))
   title('Mesh'), axis image
figure('Name','Velocity field')
   quiver(p(:,1),p(:,2),u1,u2)
   title('Velociy field'), axis image
figure('Name','Isobars')
   kpde2dcont(p,t,pr,30)
   title('Isobars'), axis image
figure('Name','Streamline')
   kpde2dcont(p,t,w,10)
   title('Streamlines'), axis image

