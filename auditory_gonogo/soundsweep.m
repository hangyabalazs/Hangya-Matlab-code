function soundsweep

octps    = 3;

wbot     = 1024;

wtop     = 1024;

wmid = sqrt(wbot*wtop);

sigmaoct = 1;

 

 

srate = 44000;

t = 0:(1/srate):2;

dt = 1/srate;

 

fulls  = zeros(1, length(t)); 

mulogw = zeros(1, length(t));

ws     = [];

if octps >= 0,  logws  = log2(wbot) - 8*sigmaoct : log2(wtop) + 4*sigmaoct;

else            logws  = log2(wbot) - 2*sigmaoct : log2(wtop) + 3*sigmaoct;

end;

 

wtrack = zeros(length(logws), length(t));

vtrack = zeros(length(logws), length(t));

 

logbot = log2(wbot); 

logtop = log2(wtop);

% vols   = [(1+tanh((logws-logbot)/sigmaoct))/2] .* ...

%         [(1+tanh((logtop-logws)/sigmaoct))/2];

 

vols = exp(-(logws-log2(wmid)).^2/(2*sigmaoct.^2));

phis   = zeros(size(logws));

 

for i=1:length(t),

   phis = phis + 2.^logws*dt;

   fulls(i)  = sum(vols.*sin(2*pi*phis));

   mulogw(i) = mean(vols.*logws);

   wtrack(:,i) = logws';

   vtrack(:,i) = vols';

   

   logws = logws + octps*dt;

   if octps >= 0  &&  logws(end) > log2(wtop) + 5.5*sigmaoct,

      logws = logws - 1;

      startw = log2(wbot) - 3*sigmaoct - 1;

      newphi = (1/(log(2)*octps))*(exp(log(2)*octps*t(i)) - 1);

      newphi = newphi * (2.^startw);

      phis  = [newphi  phis(1:end-1)];

   end;

   if octps < 0  &&  logws(1) < log2(wbot) - 2.5*sigmaoct,

      logws = logws + 1;

      phis  = [phis(2:end) 0];

   end;

   % vols   = [(1+tanh((logws-logbot)/sigmaoct))/2] .* ...

   % [(1+tanh((logtop-logws)/sigmaoct))/2];

   vols = exp(-(logws-log2(wmid)).^2/(2*sigmaoct.^2));

end;   

 

 

return;

 

for w0 = [0.5 1 2 4 8 16 32 64 128]*4,

   l2w = log2(w0) + t*octps;

   w = 2.^l2w;

   phi = cumsum(w)*dt;

   s = sin(2*pi*phi);

 

   fulls = fulls + s;

   ws = [ws ; w];

end;