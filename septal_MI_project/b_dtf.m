function b_dtf

% Import and discrimination
b_in3
b_disc
b_var2ws('all','caller')

% Sinc convolution
zint = sincconv(vdisc,length(unit));

% Estimate optimal model order
en = eeg(1:10:end)';
un = zint';
if length(en) > length(un)
    en = en(1:end-1);
elseif length(un) > length(en)
    un = un(1:end-1);
end
en = (en - mean(en)) / std(en);      % normalization
un = (un - mean(un)) / std(un);
s = [en un];
p = arorder(s,1,512);

% Calculate DTF and PDC
N = 512;    % number of frequencies
Fs = 1000;  % sampling rate
% [M,SE,R] = tfmvar(s,p,N,Fs);
[AR,RC,PE] = mvar(s,p,5);
X.A = [eye(size(AR,1)),-AR]; 
X.B = eye(size(X.A,1));
X.C = eye(size(X.A,1));
X.datatype = 'MVAR';
[S,h,PDC,COH,DTF] = main(X,'DTF',N,Fs);

% X.C = PE(:,p*size(s',1)+(1:size(s',1)))';
% [S,h,PDC,COH,DTF,DC,pCOH,dDTF,ffDTF,pCOH2,coh] = mvfreqz(X.B,X.A,X.C,N,Fs);
% 
% % Plot
% [K1,K2] = size(X.B);
% f = (1:N)/N/2*Fs;
% 
% for k=1:K1,
%     Label{k}=sprintf('#%02i',k);
% end;
% 
% figure
% for k1=1:K1;
%     for k2=1:K2;
%         subplot(K1,K2,k2+(k1-1)*K1);
%         area(f,squeeze(PDC(k1,k2,:)));
%         axis([0,max(f),0,1]);
%         if k2==1;
%             ylabel(Label{k1});
%         end;
%         if k1==1;
%             title(Label{k2});
%         end;
%     end;
% end;
% 
% figure
% for k1=1:K1;
%     for k2=1:K2;
%         subplot(K1,K2,k2+(k1-1)*K1);
%         area(f,squeeze(DTF(k1,k2,:)));
%         axis([0,max(f),0,1]);
%         if k2==1;
%             ylabel(Label{k1});
%         end;
%         if k1==1;
%             title(Label{k2});
%         end;
%     end;
% end;



% -------------------------------------------------------------------------
function zint = sincconv(vdisc,lenu)
fs = 10000;     % unit
dto = 1 / fs;
ts = zeros(1,lenu);
ts(vdisc) = 1;
du = diff(vdisc);
fdu = 1 ./ du;
fdu = [fdu 0.0001];
fcut = 100; 
fsnew = 1000;
dtnew = 1 / 1000;
fsold = 10000;
fsratio = fsnew / fsold;
told = vdisc * dto * fcut;
tnew = (1:lenu*fsratio) * dtnew * fcut;
lentold = length(told);
zint = 0;
for i = 1:lentold
    zint = zint + sinc(tnew-told(i));
end



% -------------------------------------------------------------------------
function popt = arorder(v,pmin,pmax)
% See ARFIT for a detailed help.

% Input argument check
error(nargchk(3,3,nargin))
if (pmin ~= round(pmin) | pmax ~= round(pmax))
    error('Order must be integer.');
end
if (pmax < pmin)
    error('PMAX must be greater than or equal to PMIN.')
end

% n: number of observations; m: dimension of state vectors
[n,m] = size(v);

mcor = 1;               % fit intercept vector
selector = 'sbc';	          % use SBC as order selection criterion

ne = n-pmax;               % number of block equations of size m
npmax = m*pmax+mcor;          % maximum number of parameter vectors of length m

if (ne <= npmax)
    error('Time series too short.')
end

% Compute QR factorization for model of order pmax
[R, scale] = arqr(v, pmax, mcor);

% compute approximate order selection criteria for models
% of order pmin:pmax
[sbc, fpe] = arord(R, m, mcor, ne, pmin, pmax);

% Get index iopt of order that minimizes the order selection
% criterion specified by the variable selector
[val, iopt] = min(eval(selector));

% Select order of model
popt = pmin + iopt-1; % estimated optimum order
np = m*popt + mcor; % number of parameter vectors of length m



% -------------------------------------------------------------------------
function [S,h,PDC,COH,DTF] = main(X,Mode,N,Fs)
% See PLOTA for a detailed help.

% Initialize
[K1,K2] = size(X.A);
p = K2/K1-1;
[K1,K2] = size(X.B);
q = K2/K1-1;
f = (1:N)/N/2*Fs;
z = i*2*pi/Fs;

h = zeros(K1,K1,N);
S = zeros(K1,K1,N);
DTF = zeros(K1,K1,N);
COH = zeros(K1,K1,N);
PDC = zeros(K1,K1,N);
PDCF = zeros(K1,K1,N);
invC = inv(X.C);
tmp1 = zeros(1,K1);
tmp2 = zeros(1,K1);

% Calculate PDC
for n = 1:N
    atmp = zeros(K1,K1);
    for k = 1:p+1,
        atmp = atmp + X.A(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
    end
    btmp = zeros(K1,K2);
    for k = 1:q+1,
        btmp = btmp + X.B(:,k*K1+(1-K1:0))*exp(z*(k-1)*f(n));
    end
    h(:,:,n) = atmp\btmp;
    S(:,:,n) = h(:,:,n)*X.C*h(:,:,n)';

    for k1 = 1:K1,
        tmp = squeeze(atmp(:,k1));
        tmp1(k1) = sqrt(tmp'*tmp);
        tmp2(k1) = sqrt(tmp'*invC*tmp);
    end

    PDCF(:,:,n) = abs(atmp)./tmp2(ones(1,K1),:)	;
    PDC(:,:,n)  = abs(atmp)./tmp1(ones(1,K1),:);
end

% Calculate DTF
for k1=1:K1;
    DEN=sqrt(sum(abs(h(k1,:,:)).^2,2));
    for k2=1:K2;
        COH(k1,k2,:) = abs(S(k1,k2,:))./sqrt(abs(S(k1,k1,:).*S(k2,k2,:)));
        DTF(k1,k2,:) = abs(h(k1,k2,:))./DEN;
    end
end

% Plot
for k=1:K1,
    Label{k}=sprintf('#%02i',k);
end;

figure
for k1=1:K1;
    for k2=1:K2;
        subplot(K1,K2,k2+(k1-1)*K1);
        area(f,squeeze(PDC(k1,k2,:)));
        axis([0,max(f),0,1]);
        if k2==1;
            ylabel(Label{k1});
        end;
        if k1==1;
            title(Label{k2});
        end;
    end;
end;

figure
for k1=1:K1;
    for k2=1:K2;
        subplot(K1,K2,k2+(k1-1)*K1);
        area(f,squeeze(DTF(k1,k2,:)));
        axis([0,max(f),0,1]);
        if k2==1;
            ylabel(Label{k1});
        end;
        if k1==1;
            title(Label{k2});
        end;
    end;
end;



% -------------------------------------------------------------------------
function [M,SE,R] = tfmvar(s,MOP,f,Fs)
% See TFMVAR for a detailed help.

FLAG.SEM = 1; 

R.datatype = 'TF-MVAR';
R.N = length(TRIG);
R.T = (T(1,:)+T(2,:))/(2*Fs); 
R.F = f; 
R.MOP = MOP; 
R.nan_ratio = zeros(1,size(T,2));


% univariate
[N,m] = size(s);

sz = [m,length(f),size(T,2)];
R.M.S1 = zeros(sz);
R.M.logS1   = zeros(sz);

sz(2)  = MOP;
R.M.AR1 = zeros(sz);
sz(2)  = MOP+1;
R.M.C1  = zeros(sz);

if FLAG.SEM,
        R.SE1 = R.M;
end;

sz    = [m,length(f),length(TRIG)];
S1     = zeros(sz);
sz(2) = MOP;
A1 = zeros(sz);
sz(2) = MOP+1;
C1  = zeros(sz);

% multivariate 
[N,m] = size(s);

sz = [m,m,length(f),size(T,2)];
R.M.S = zeros(sz);
R.M.h = zeros(sz);
R.M.phaseS = zeros(sz);
R.M.phaseh = zeros(sz);
R.M.logS   = zeros(sz);
R.M.logh   = zeros(sz);
R.M.COH    = zeros(sz);
R.M.coh    = zeros(sz);
R.M.PDC    = zeros(sz);
R.M.DTF    = zeros(sz);
R.M.pCOH   = zeros(sz);
R.M.dDTF   = zeros(sz);
R.M.ffDTF  = zeros(sz);
R.M.pCOH2  = zeros(sz);
sz(3)  = 1;
R.M.DC = zeros(sz);
R.M.C  = zeros(sz);
sz(2)  = m*MOP;
R.M.AR = zeros(sz);

R.SE = R.M;

sz    = [m,m,length(f),length(TRIG)];
S     = zeros(sz);
h     = zeros(sz);
COH   = zeros(sz);
coh   = zeros(sz);
PDC   = zeros(sz);
DTF   = zeros(sz);
pCOH  = zeros(sz);
dDTF  = zeros(sz);
ffDTF = zeros(sz);
pCOH2 = zeros(sz);
sz(3) = 1;
C  = zeros(sz);
DC = zeros(sz);
sz(2) = m*MOP;
AR = zeros(sz);

tic,
for k1 = 1:size(T,2),
        [S0,sz0] = trigg(s,TRIG,T(1,k1),T(2,k1),1+MOP);

        R.nan_ratio(k1) = mean(isnan(S0(:)));
        % univariate
        [AR1,RC1,PE1] = durlev(acovf(S0,MOP));
        
        for k=1:m;
                [h1,f] = freqz(sqrt(PE1(k,MOP+1)/(Fs*2*pi)),ar2poly(AR1(k,:)),f,Fs);
                H1(k,:)= h1(:)'; %F(:,k)=f(:);
        end;
        
        %[S(:,:,:,k2),  h(:,:,:,k2), PDC(:,:,:,k2), COH(:,:,:,k2), DTF(:,:,:,k2), DC(:,:,1,k2), pCOH(:,:,:,k2), dDTF(:,:,:,k2), ffDTF(:,:,:,k2), pCOH2(:,:,:,k2),coh(:,:,:,k2)] = mvfreqz(X.B,X.A,X.C,f,Fs);
        
        R.M.S1(:,:,k1) = H1; 
        R.M.logS1(:,:,k1) = log(abs(H1)); 
        R.M.AR1(:,:,k1) = AR1; 
        R.M.C1(:,:,k1)  = PE1;
        
        % multivariate 
        [A,RC,PE] = mvar(S0',MOP,5);

        X.A = [eye(size(S0,1)),-A];

        X.B = [eye(size(S0,1))];

        X.C = PE(:,MOP*size(S0,1)+(1:size(S0,1)))';

        X.datatype = 'MVAR';

        
        [S,  h, PDC, COH, DTF, DC, pCOH, dDTF, ffDTF, pCOH2,coh] = mvfreqz(X.B,X.A,X.C,f,Fs);

        
        R.M.phaseS(:,:,:,k1) = angle(S);
        R.M.phaseh(:,:,:,k1) = angle(h);
        R.M.S(:,:,:,k1)      = S;
        R.M.h(:,:,:,k1)      = h;
        R.M.logS(:,:,:,k1)   = log(abs(S));
        R.M.logh(:,:,:,k1)   = log(abs(h));
        R.M.PDC(:,:,:,k1)    = PDC;
        R.M.COH(:,:,:,k1)    = COH;
        R.M.coh(:,:,:,k1)    = coh;
        R.M.DTF(:,:,:,k1)    = DTF;
        R.M.pCOH(:,:,:,k1)   = pCOH;
        R.M.dDTF(:,:,:,k1)   = dDTF;
        R.M.ffDTF(:,:,:,k1)  = ffDTF;
        R.M.pCOH2(:,:,:,k1)  = pCOH2;
        R.M.AR(:,:,1,k1)     = A;
        R.M.DC(:,:,1,k1)     = DC;
        R.M.C(:,:,1,k1)      = X.C;
        
        if FLAG.SEM,
                for k2 = 1:length(TRIG),
                        sel = [1:k2-1,k2+1:length(TRIG)];
                        [S0,sz0] = trigg(s, TRIG(sel), T(1,k1), T(2,k1), 1+MOP);

                        %fprintf(2,'\nExtract epochs done\nNumber of Trials: %i : %i\n',length(t0),length(t1));

                        

                        % univariate
                        [AR1,RC1,PE1] = durlev(acovf(S0,MOP));
                        
                        for k=1:m;
                                [h1,f] = freqz(sqrt(PE1(k,MOP+1)/(Fs*2*pi)),ar2poly(AR1(k,:)),f,Fs);
                                H1(k,:)= h1(:)'; %F(:,k)=f(:);
                        end;
                        
                        %[S(:,:,:,k2),  h(:,:,:,k2), PDC(:,:,:,k2), COH(:,:,:,k2), DTF(:,:,:,k2), DC(:,:,1,k2), pCOH(:,:,:,k2), dDTF(:,:,:,k2), ffDTF(:,:,:,k2), pCOH2(:,:,:,k2),coh(:,:,:,k2)] = mvfreqz(X.B,X.A,X.C,f,Fs);
                        
                        S1(:,:,k2) = H1; 
                        A1(:,:,k2) = AR1; 
                        C1(:,:,k2)  = PE1;
                        
                        % multivariate 
                        [A,RC,PE] = mvar(S0',MOP,5);

                        X.A = [eye(size(S0,1)),-A];

                        X.B = [eye(size(S0,1))];

                        X.C = PE(:,MOP*size(S0,1)+(1:size(S0,1)))';

                        X.datatype = 'MVAR';

                        
                        [S(:,:,:,k2),  h(:,:,:,k2), PDC(:,:,:,k2), COH(:,:,:,k2), DTF(:,:,:,k2), DC(:,:,1,k2), pCOH(:,:,:,k2), dDTF(:,:,:,k2), ffDTF(:,:,:,k2), pCOH2(:,:,:,k2),coh(:,:,:,k2)] = mvfreqz(X.B,X.A,X.C,f,Fs);

                        
                        AR(:,:,1,k2) = A; 
                        C(:,:,1,k2)  = X.C;
                        
                        fprintf(2,'%.0f\t%.0f\t%.1f\n',[k1,k2,toc]);tic;
                end;
                
                % univariate
                [R.SE.S1(:,:,k1),   R.M.S1(:,:,k1)  ] = sem(abs(S1),3);
                [R.SE.logS1(:,:,k1),R.M.logS1(:,:,k1)] = sem(log(abs(S1)),3);
                [R.SE.AR1(:,:,k1),  R.M.AR1(:,:,k1) ] = sem(A1,3);
                [R.SE.C1(:,:,k1),   R.M.C1(:,:,k1)  ] = sem(C1,3);
                
                
                % multivariate 
                [R.SE.phaseS(:,:,:,k1), R.M.phaseS(:,:,:,k1)] = sem(angle(S),4);
                [R.SE.phaseh(:,:,:,k1), R.M.phaseh(:,:,:,k1)] = sem(angle(h),4);
                [R.SE.S(:,:,:,k1),      R.M.S(:,:,:,k1)     ] = sem(S,4);
                [R.SE.h(:,:,:,k1),      R.M.h(:,:,:,k1)     ] = sem(h,4);
                [R.SE.logS(:,:,:,k1),   R.M.logS(:,:,:,k1)  ] = sem(log(abs(S)),4);
                [R.SE.logh(:,:,:,k1),   R.M.logh(:,:,:,k1)  ] = sem(log(abs(h)),4);
                [R.SE.PDC(:,:,:,k1),    R.M.PDC(:,:,:,k1)   ] = sem(PDC,4);
                [R.SE.COH(:,:,:,k1),    R.M.COH(:,:,:,k1)   ] = sem(COH,4);
                [R.SE.coh(:,:,:,k1),    R.M.coh(:,:,:,k1)   ] = sem(coh,4);
                [R.SE.DTF(:,:,:,k1),    R.M.DTF(:,:,:,k1)   ] = sem(DTF,4);
                [R.SE.pCOH(:,:,:,k1),   R.M.pCOH(:,:,:,k1)  ] = sem(pCOH,4);
                [R.SE.dDTF(:,:,:,k1),   R.M.dDTF(:,:,:,k1)  ] = sem(dDTF,4);
                [R.SE.ffDTF(:,:,:,k1),  R.M.ffDTF(:,:,:,k1) ] = sem(ffDTF,4);
                [R.SE.pCOH2(:,:,:,k1),  R.M.pCOH2(:,:,:,k1) ] = sem(pCOH2,4);
                [R.SE.AR(:,:,1,k1),     R.M.AR(:,:,1,k1)    ] = sem(AR,4);
                [R.SE.DC(:,:,1,k1),     R.M.DC(:,:,1,k1)    ] = sem(DC,4);
                [R.SE.C(:,:,1,k1),      R.M.C(:,:,1,k1)     ] = sem(C,4);
                
        end; 
        %% In order to obtaine the standard error of the mean, 
        %% this SE must be multiplied by N-1, with N = length(TRIG)
end;

if any(R.nan_ratio),
        fprintf(2,'Warning TFMVAR: up to %5.2f%% are missing values (NaN).\n',max(R.nan_ratio)*100);
        fprintf(2,'This may cause underestimating the standard error.\n');
end;	

if nargout > 1, 
        M = R.M;
        R = rmfield(R,'M');
        SE= R.SE; 
        R = rmfield(R,'SE');
else 
        M = R;
end;