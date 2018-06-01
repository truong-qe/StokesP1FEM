%
% 3D Stokes Problem with P1-Bubble/P1 FEM
%
clear all, close all
%
alpha=100;
nu=0.01;
%
% choose method 1: LDL factorization, 2: Uzawa Conjugate Gradient
meth=2;   

% load mesh
load cube4913 p t ibc ibcp 
np=size(p,1); nt=size(t,1);  
ibcp=1; 
fprintf('Nodes=%5d  Triangles=%5d \n',np,nt)

x=p(:,1); y=p(:,2); z=p(:,3);
u1e=y; u2e=z; u3e=x; pre=x+y+z;
f1=ones(np,1)+alpha*u1e; f2=ones(np,1)+alpha*u2e; f3=ones(np,1)+alpha*u3e;
u1bc=u1e(ibc); u2bc=u2e(ibc); u3bc=u3e(ibc); 
 
%
% Assemble matrices and right-hand side
if (meth == 1)
   A=kstok3dp1bmat(p,t,nu,alpha); 
   b=kstok3dp1brhs(p,t,f1,f2,f3,nu,alpha);
elseif (meth == 2)
    [K,B,C]=kstok3dp1bmat(p,t,nu,alpha); 
    [b,bp]=kstok3dp1brhs(p,t,f1,f2,f3,nu,alpha);
    A1=K; A2=A1; A3=A1;
    
    % Preconditioner
    K0=kpde3dstf(p,t,1); M0=kpde3dmss(p,t,1);
    DM0=spdiags(M0,0); DM0=spdiags(DM0,0,np,np);
    vol=kpde3dgphi(p,t);
    hx=1/16; Rnh=alpha*hx*hx/nu;
    Kh=K0+alpha*C; KM=nu*K0; 
    if (Rnh <= 1)
       KM=KM+alpha*DM0;
    else
       KM=KM+alpha*speye(np);
    end
end

tic
%
% Incorporate the boundary conditions
PenBC=10^15;
if (meth == 1)
   ibc1=ibc; ibc2=np+ibc; ibc3=2*np+ibc; ibcp=3*np+1;
   A(ibc1,ibc1)=A(ibc1,ibc1)+PenBC*speye(length(ibc1)); 
   A(ibc2,ibc2)=A(ibc2,ibc2)+PenBC*speye(length(ibc2)); 
   A(ibc3,ibc3)=A(ibc3,ibc3)+PenBC*speye(length(ibc3));
   A(ibcp,ibcp)=A(ibcp,ibcp)+PenBC*speye(length(ibcp));
   b(ibc1)=PenBC*u1bc; b(ibc2)=PenBC*u2bc;  b(ibc3)=PenBC*u3bc; b(ibcp)=0;
elseif (meth == 2) 
   ibc1=ibc; ibc2=ibc; ibc3=ibc; ibcp=1;
   A1(ibc1,ibc1)=A1(ibc1,ibc1)+PenBC*speye(length(ibc1)); 
   A2(ibc2,ibc2)=A2(ibc2,ibc2)+PenBC*speye(length(ibc2)); 
   A3(ibc3,ibc3)=A3(ibc3,ibc3)+PenBC*speye(length(ibc3));
   Kh(ibcp,ibcp)=Kh(ibcp,ibcp)+PenBC*speye(length(ibcp));
   b(ibc1)=PenBC*u1bc; b(np+ibc2)=PenBC*u2bc; b(2*np+ibc3)=PenBC*u3bc; bp(ibcp)=0;
end       

%
% Compute the solution (LDL' factorization)
if (meth == 1)
    [L,D,s]=ldl(A,'vector');
    u=zeros(4*np,1); u(s)=L'\(D\(L\b(s)));
    u1=u(1:np); u2=u(np+1:2*np); u3=u(2*np+1:3*np); pr=u(3*np+1:end);
    fprintf('LDL factorization CPU=%12.6f \n',toc)
elseif (meth == 2)
    [R1,p1,s1]=chol(A1,'lower','vector');
    [R2,p2,s2]=chol(A2,'lower','vector');
    [R3,p3,s3]=chol(A3,'lower','vector');
    [Rh,pp,sh]=chol(Kh,'lower','vector');
    
    % Data structure
    R={R1,R2 R3}; s={s1 s2 s3}; ibcd={ibc1 ibc2 ibc3};  
    b={b(1:np) b(np+1:2*np) b(2*np+1:3*np)};
    [u,pr,iter]=kstokcg(R,Rh,B,C,KM,s,sh,b,bp,ibcd,ibcp,10^-6);
    u1=u(1:np); u2=u(np+1:2*np); u3=u(2*np+1:3*np);
    fprintf('Uzawa Conjugate Gradient : iter=%5d CPU=%12.6f\n',iter,toc)
end

 
%
% Outputs
kstokshow(p,t,u1,u2,u3,pr,ibc) 
