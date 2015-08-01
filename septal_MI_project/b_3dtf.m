function b_3dtf

% Import and discrimination
b_in3
b_disc
b_var2ws('all','caller')

% Initialize global LIN
global LIN1
global LIN2
LIN1 = [];
LIN2 = [];

% Main
seglen = 20000;
maxi = floor(length(eeg)/seglen);
PDC = zeros(2,2,512,maxi);      % for N = 512 !
for i = 1:maxi
    ind1 = (i - 1) * seglen + 1;
    ind2 = ind1 + seglen -1;
    eeg2 = eeg(ind1:ind2);
    unit2 = unit(ind1:ind2);
    vdisc2 = vdisc(find(vdisc>=ind1&vdisc<=ind2)) - ind1;
    
    if ~isempty(vdisc2)
        PDC(:,:,:,i) = dircoh(eeg2,unit2,vdisc2);  % partial directed coherence
    else
        PDC(:,:,:,i) = 0;
    end
end
% mPDC = squeeze(max(PDC,[],3));

% Plot
K1 = 2;
K2 = 2;
for k = 1:K1
    Label{k} = sprintf('#%02i',k);
end

% figure
% for k1 = 1:K1
%     for k2 = 1:K2
%         subplot(K1,K2,k2+(k1-1)*K1);
%         imagesc(squeeze(PDC(k1,k2,:)));
% %         axis([0,maxi*seglen,0,1]);
%         if k2 == 1;
%             ylabel(Label{k1});
%         end
%         if k1 == 1;
%             title(Label{k2});
%         end
%     end
% end

% Plot LIN
global LIN1
global LIN2
figure
plot(LIN1,'r')  % unit->eeg
hold on
plot(LIN2,'b')  % eeg->unit
% plot(LIN1-LIN2,'g')



% -------------------------------------------------------------------------
function PDC = dircoh(eeg,unit,vdisc)

% Blackman-Harris convolution
if vdisc(1) == 0
    vdisc = vdisc(2:end);
end
zint = bhconv(vdisc,length(unit));
zint = zint(1:10:end);

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
% en = en ./ max(en);      % normalization
% un = un ./ max(un);
s = [en un];
p = arorder(s,1,512);

% Calculate DTF and PDC
N = 512;    % number of frequencies
Fs = 1000;  % sampling rate
[AR,RC,PE] = mvar(s,p,5);
X.A = [eye(size(AR,1)),-AR]; 
X.B = eye(size(X.A,1));
X.C = eye(size(X.A,1));
X.datatype = 'MVAR';
[S,h,PDC,COH,DTF] = main(X,'DTF',N,Fs);



% -------------------------------------------------------------------------
function zint = bhconv(vdisc,lenu)
ipunit = zeros(1,lenu);
ipunit(vdisc) = 1;
wbh = blackmanharris(1000);
wipunit = conv(ipunit,wbh);
zint = wipunit(1:lenu);



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

% Ask global LIN
global LIN1
global LIN2

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
%         axis([0,100,0,1]);
        sqPDC = squeeze(PDC(k1,k2,:));
        tt1 = num2str(sum(sqPDC));
        tt2 = num2str(max(sqPDC));
        text(350,0.9,tt1)
        text(350,0.6,tt2)
        if k2==1;
            ylabel(Label{k1});
        end;
        if k1==1;
            title(Label{k2});
        end;
        if k1 == 1 & k2 == 2
            LIN1(end+1) = sum(sqPDC(3:6));
        end
        if k1 == 2 & k2 == 1
            LIN2(end+1) = sum(sqPDC(3:6));
        end
    end;
end;

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