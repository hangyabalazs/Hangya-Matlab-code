function [x1,x2,hx,hy,jh,N,Hx,Hy,Hxy,Hycx,Hxcy,Ixy,Uxy,Uyx,Ux2y,m,n] = entr3(tfv1,tfv2,y1,l1,y2,l2)

 n1=length(y1);
 n2=length(y2);

 % histogram estimation
 
    miny1 = min(min(tfv1));               %min/max of data: subject to be changed
    maxy1 = max(max(tfv1));
    binwidth1 = (maxy1 - miny1) ./ l1;
    xx1 = miny1 + binwidth1*(0:l1);
    xx1(length(xx1)) = maxy1;
    x1 = xx1(1:length(xx1)-1) + binwidth1/2;


    miny2 = min(min(tfv2));               %min/max of data: subject to be changed
    maxy2 = max(max(tfv2));
    binwidth2 = (maxy2 - miny2) ./ l2;
    xx2 = miny2 + binwidth2*(0:l2);
    xx2(length(xx2)) = maxy2;
    x2 = xx2(1:length(xx2)-1) + binwidth2/2;

nbin1 = length(xx1);
nbin2 = length(xx2);
jh= zeros(nbin1,nbin2);


for i=2:nbin1
   if i==2
     ix1=find((y1<=xx1(i)));
   else
     ix1=find((y1<=xx1(i))&(y1>xx1(i-1)));
   end
   for k=2:nbin2
      if k==2
        ix2=find((y2<=xx2(k)));
      else
        ix2=find((y2<=xx2(k))&(y2>xx2(k-1)));
     end
     ix=[];
     if ((isempty(ix1))|(isempty(ix2)))==0		
      lenix1=length(ix1);
      lenix2=length(ix2);
      if lenix1>lenix2
         lenix=lenix1;
         ix121=ix1;
         ix122=ix2;
      else
         lenix=lenix2;
         ix121=ix2;
         ix122=ix1;
      end
      
      for m=1:lenix
         if (isempty(ix121))|(isempty(ix122))==1 
            disp('NOW!!')
            keyboard
         end
         vtemp=find(ix121(m)==ix122);
         if isempty(vtemp)~=1
            ix=[ix vtemp ];
        end
       end
    end
    % keyboard
   jh(i,k)=length(ix);
   end
end

% calculation of entropies & uncertainity coefficients

[m,n]=size(jh); 
N=sum(sum(jh));


hxx=sum(jh); hyy=sum(jh');

N1=sum(hxx);
N2=sum(hyy);
hx=hxx/N1;
hy=hyy/N2;

Hx=0;
for i=1:m
   if hx(i)~=0
      a=hxx(i)/N;
      Hx=Hx-(a*log(a));
   end
end

Hy=0;
for k=1:n
   if hy(k)~=0
      a=hyy(k)/N;
      Hy=Hy-(a*log(a));
   end
end

Hxy=0;
for i=1:m
   for k=1:n 
      if jh(i,k)~=0
         a=jh(i,k)/N;
         Hxy=Hxy-a*log(a);
      end
   end
end


Hycx=Hxy-Hx;
Hxcy=Hxy-Hy;

tiny=1.0e-30;
Uyx=(Hy-Hycx)/(Hy+tiny);
Uxy=(Hx-Hxcy)/(Hx+tiny);
Ux2y=2*((Hy+Hx-Hxy)/(Hx+Hy+tiny));

Ixy=Hx+Hy-Hxy;

entropies=[Hx,Hy,Hxy,Hycx,Hxcy,Ixy,Uxy,Uyx,Ux2y];
[k1,k2]=size(entropies);
if k1<k2
   entropies=entropies';
end





