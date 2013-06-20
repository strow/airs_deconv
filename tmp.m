
>> Binv = pinv(full(B1.Cmat(A1.sind,:))) ;
>> B = full(B1.Cmat(A1.sind,:)) ;         

>> whos B Binv   
  Name          Size                   Bytes  Class     Attributes

  B          2314x42860            793424320  double              
  Binv      42860x2314             793424320  double              

>> plot(vc1d, 50*B(iy, :), vc1d, Binv(:,iy), vc1d, 100*Bsum-40)
>> plot(vc1d, 50*B(iy, :), vc1d, Binv(:,iy), vc1d, 100*Bsum-10)
>> plot(vc1d, 50*B(iy, :), vc1d, Binv(:,iy), vc1d, 100*Bsum-11)
>> grid
>> zoom on
>> zz = sum(Binv);
>> whos zz
  Name      Size              Bytes  Class     Attributes

  zz        1x2314            18512  double              

>> zz = sum(Binv,2);
>> plot(zz)
>> whos zz          
  Name          Size             Bytes  Class     Attributes

  zz        42860x1             342880  double              

