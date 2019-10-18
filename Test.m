clear all;close all
nr=326;nc=351;
load XYZ;
load XI;
load YI;

%XI=XI*xscale;YI=YI*yscale;
map=zeros(nr,nc);
n=length(XYZ(:,1));
xx=XYZ(:,1);yy=XYZ(:,2);zz=XYZ(:,3);
[xx1,yy1,zz1] = mincurvi(xx,yy,zz,XI,YI);
xyz1=zeros(nr,nc);
k=0;
for i=1:nr
    for j=1:nc
        k=k+1;
        xyz1(i,j)=zz1(k);
    end
end
figure;contour(xyz1);

