function [u,p,iter] = kstokcg(R,Rp,B,C,H,s,sp,b,bp,ibc,ibcp,tol)
% 
%
itmax=300;
np=length(bp); nd=length(b); n=nd*np;
p=zeros(np,1); gk=zeros(np,1); g=zeros(np,1); z=zeros(np,1);
% Initialization
v=cell(nd,1);  [v{:}]=deal(zeros(np,1));
for i=1:nd, v{i}(s{i})=R{i}'\(R{i}\b{i}(s{i})); end
u=cell2mat(v);
r=B*u+C*p+bp; bp=H*r; bp(ibcp)=0; g(sp)=Rp'\(Rp\bp(sp));  
d=g; 
gg=g'*r; gg0=max(1,gg); 

iter=0; err=1;
%fprintf('iter=%5d  |g|=%15.8e  |r|=%15.8e |g|_2=%15.8e \n',iter,sqrt(gg),sqrt(r'*r),sqrt(g'*g)) 
%pause
while (err>tol && iter<itmax)
    iter=iter+1;
    % Sensitivity
    bb=B'*d;
    for i=1:nd
        b{i}=bb(np*(i-1)+1:i*np);
        b{i}(ibc{i})=0; v{i}(s{i})=R{i}'\(R{i}\b{i}(s{i}));
    end
    w=cell2mat(v);
    % Step size
    rk=B*w+C*d; rk(ibcp)=0; 
    tk=(r'*d)/(rk'*d); 
    % Update
    p=p-tk*d; u=u-tk*w; r=r-tk*rk;
    % New gradient
    bp=H*rk;  bp(ibcp)=0;  gk(sp)=Rp'\(Rp\bp(sp)); 
    g=g-tk*gk;
    % New direction
    gg1=g'*r;  beta=gg1/gg; 
    d=g+beta*d;
    % Stopping criterion
    err=sqrt(gg1/gg0); 
    %fprintf('iter=%5d  tk=%15.8e |r|=%15.8e  err=%15.8e \n',iter,tk,sqrt(r'*r),err) 
    gg=gg1;
end

