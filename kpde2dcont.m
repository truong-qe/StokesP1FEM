function []=kpde2dcont(p,t,u,n,umin,umax)
%KPDE2DCONT Draws n level curves of the approximate solution u 
%--------------------------------------------------------------------
% kpde2dcont(p,t,u,n)
% Input:
%      p : Node coodinates, np*2
%      t : Triangle vertices, nt*3
%      u : Approximate solution, np*1
%      n : Number of level curves
%--------------------------------------------------------------------
% (c) J. Koko, ISIMA 03/2006, koko@isima.fr
%--------------------------------------------------------------------
% 1. Level curve values
if (nargin == 4)
   umin=min(u);  umax=max(u);
end
du=(umax-umin)/(n-1);
ucn=[umin:du:umax]'; ncn=length(ucn);

% 2. Boundary edges
a=[t(:,[1 2]); t(:,[2 3]); t(:,[1 3])]; a=sort(a,2);
[aa,ii,jj]=unique(a,'rows');
ca=histc(jj,1:max(jj)); i1=find(ca==1); abord=a(ii(i1),:);
xb=[p(abord(:,1),1) p(abord(:,2),1)];
yb=[p(abord(:,1),2) p(abord(:,2),2)];

% 3. Loop over triangles
ns=0;
nt=size(t,1);
for ih=1:nt
   it=t(ih,1:3)'; ut=u(it); pt=p(it,:);
   for il=1:ncn
      uc=ucn(il);
      j1=find(ut<=uc); j2=find(ut>=uc);
      
      % level curve does not intersect triangle ih 
      if (length(j1)==0 | length(j2)==0)
         continue
      end
      
      % level curve intersects through triangle ih
      if (length(j1)==1)
         i1=j1; i2=j2(1); i3=j2(2);
      else
         i1=j2; i2=j1(1); i3=j1(2);
      end
      a1=pt(i2,:)-pt(i1,:); a2=pt(i3,:)-pt(i1,:);
      
      % compute intersction edges/curve points
      z1=sqrt(sum(a1.^2,2)); z2=sqrt(sum(a2.^2,2));
      t1=(uc-ut(i1)+eps)/(ut(i2)-ut(i1)+eps);
      t2=(uc-ut(i1)+eps)/(ut(i3)-ut(i1)+eps);
      p1=pt(i1,:)+t1*a1; p2=pt(i1,:)+t2*a2;
      
      % level curve in triangle ih
      ns=ns+1;
      xc(1,ns)=p1(1); xc(2,ns)=p2(1);
      yc(1,ns)=p1(2); yc(2,ns)=p2(2);
   end      
end

% 4. Draw level curves (blue) and the boundary edges (black)
% plot(xb',yb','k-',xc,yc,'b-','LineWidth',1.75)
plot(xb',yb','k-',xc,yc,'b-')

