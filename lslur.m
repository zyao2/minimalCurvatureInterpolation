function B=lslur(N,B,A)
for I = 1:N
  IK =(I*(I-1)/2);
  for J = 1:I-1
     JK = (J*(J-1)/2);
     IKJ = IK+J;
     sum=0;
     for K = 1:J-1
       IKK = IK+K;
       JKK = JK+K;
       sum=sum+A(IK+K)*A(JK+K);
     end  
      A(IKJ) = A(IKJ)-sum;
   end   
   IKI = IK+I;
   sum=0;
   for K = 1: I-1
      IKK = IK+K;
      Y = (A(IKK));     
      Z = Y / A(K*(K+1)/2);       
      A(IKK) = Z;
      sum=sum+ Y*Z;
   end
   A(IKI) =A(IKI)-sum;
   if( A(IKI) == 0.)
       PRINT *,'---ERR, SINGULAR MATRIX',IK,I,IK+I
       stop
   end
end

for I = 1:N
   IM1 = I - 1;
   IK = I*IM1/2;

   sum=0;
   if(I>1)
      for K = 1:IM1
         IKK = IK+K;
         sum=sum+ A(IKK)*B(K);      
      end
    end     
    B(I) =  B(I)-sum;
 end
 
 for I = 1:N
    B(I) = B(I) / A(I*(I+1)/2);
 end     
 for I = N:-1:1
    IM1 = I - 1;
    IK = I*IM1/2;
    X =(B(I));
    if(I>1)
      for K = 1:IM1
         IKK = IK+K; 
         B(K) = B(K) - X*A(IKK);
      end
  end  
end

