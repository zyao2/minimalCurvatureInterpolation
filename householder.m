function F=householder(n,B,A)
%A(n*(n+1)/2,1);
% for i=1:n
%    for j=1:i-1
%       A((i-1)*i/2+j)::A(i,j)
C=zeros(n,1);
nn=n*(n+1)/2;
n2=n-2;
for k=1:n2
    aa=0.;
	k1=k+1;
	for i=k1:n
        aa=aa+A((i-1)*i/2+k)^2;
    end
	bb=sqrt(aa);
	kk=k*(k+3)/2;
	h1=aa+abs(A(kk))*bb;
    C(k+1)=-A(kk)/abs(A(kk))*bb;
	A(kk)=A(kk)-C(k+1);
	for i=k1:n
        F(i)=0.;
	    for j=k1:n
            if(i>=j)
                ii=i;
            else
                ii=j;
            end
            jj=i+j-ii;
            F(i)=F(i)+A((ii-1)*ii/2+jj)*A((j-1)*j/2+k);
        end
   	    F(i)=F(i)/h1;
    end
	bb=0.;
	aa=0.;
	for i=k1:n
	    kk=(i-1)*i/2+k;
	    bb=bb+A(kk)*F(i);
   	    aa=aa+A(kk)*B(i);
    end
	bb=bb/(2.0*h1);
	aa=aa/h1;
	for i=k1:n
	    kk=(i-1)*i/2+k;
	    F(i)=F(i)-bb*A(kk);
  	    B(i)=B(i)-aa*A(kk);
    end
	for i=k1:n
        for j=k1:i
            kk=(i-1)*i/2 ;
            A(kk+j)=A(kk+j)-A(kk+k)*F(j)-A((j-1)*j/2+k)*F(i);
        end
    end
end

C(n)=A(nn-1);
if(abs(A(nn))<=1.E-10)
    for k=2:-1:1
        m=n-k;
        h1=B(m);
        aa=C(m);
        tt=A(m*(m+1)/2);
        for i=m:-1:3
            h1=h1-tt*B(i-1)/C(i);
            bb=aa;
            aa=-tt*C(i-1)/C(i);
            tt=bb-tt*A((i-1)*i/2)/C(i);
        end
        U(k)=(A(1)*tt/C(2)-aa)/C(m+1);
        V(k)=(h1-B(1)*tt/C(2))/C(m+1);
    end
    F(1)=(B(n)-C(n)*V(2)-A(nn)*V(1))/(C(n)*U(2)+A(nn)*U(1));
    F(2)=(B(1)-A(1)*F(1))/C(2);
    for i=2:n-1
        F(i+1)=(B(i)-C(i)*F(i-1)-A(i*(i+1)/2)*F(i))/C(i+1);
    end
else
    B(n)=B(n)/A(nn);
    A(nn)=-C(n)/A(nn);
    for i=2:n
        m=n+1-i;
        m1=m+1;
        A(m*m1/2)=C(m1)*A(m1*(m1+1)/2)+A(m*m1/2);
        B(m)=(B(m)-C(m1)*B(m1))/A(m*m1/2);
        A(m*m1/2)=-C(m)/A(m*m1/2);
    end
    F(1)=B(1);
	for i=n-1:-1:1
        m1=n+1-i;
        F(m1)=A(m1*(m1+1)/2)*F(n-i)+B(m1);
    end
end

for k=n-2:-1:1
    k1=k+1;
    bb=0.;
    h1=abs((C(k+1)*A(k*(k+3)/2)));
    for i=k1:n
        bb=bb+A((i-1)*i/2+k)*F(i);
    end
    bb=bb/h1;
	for i=k1:n
        F(i)=F(i)-bb*A((i-1)*i/2+k);
    end
end


