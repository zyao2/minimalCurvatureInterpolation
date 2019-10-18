function [xi,yi,zi] = mincurvi(xb,yb,zb,xi,yi)

is_scale = 1;     % Scale variables

damp=0.005;%if householder damp=0.;

 % Scales ...............
if is_scale
  xsc = max(xb)-min(xb);
  ysc = max(yb)-min(yb);
  xsc=1/xsc;
  ysc=1/ysc;
  xb = xb*xsc; xi = xi*xsc;
  yb = yb*ysc; yi = yi*ysc;
end


 % Calculate quadtree divison into blocks and
 % related objects


zb = zb(:);


 % Remove mean and trend ................
[zb,r0,g] = detrend2(xb,yb,zb);

zsc = max(zb)-min(zb);

 % Quick exit if residual z is 0
if zsc==0
  % Put mean and slope back
  zi = g(1)+(xi-r0(1))*g(2)+(yi-r0(2))*g(3);
  if nargout<=1, xi = zi; end
  return
end
%zb = zb/zsc;

 % Initialize output and weights
zi = zeros(size(xi));
nb=length(xb);

A=zeros(nb*(nb+1)/2,1);
for i=1:nb
    A(i*(i+1)/2)=damp;
	for j=1:i-1
        aa=xb(i)-xb(j);
        bb=yb(i)-yb(j);
        rr=aa^2+bb^2;
        if(rr>0)
            rr=rr*(log(rr)-2);
        end
        A((i-1)*i/2+j)=rr;
    end
end

v=lslur(nb,zb,A);
%v=householder(nb,zb,A);

ni=length(xi);


      
for k1=1:ni
    aa=0;
    for k2=1:nb
        Dx = xi(k1)-xb(k2);
        Dy = yi(k1)-yb(k2);
        rr = Dx.^2+Dy.^2;
        if(rr>0) 
            rr = rr.*(log(rr)-2);
        end
        aa=aa+rr*v(k2);
    end
    zi(k1) = aa;
    zi(k1) = zi(k1)+g(1)+(xi(k1)-r0(1))*g(2)+(yi(k1)-r0(2))*g(3);
end


