% from Hyun

%% load

load('C:\Balazs\_data\VIP\freq_tuning_data_inh.mat')   % inhibited group
load('C:\Balazs\_data\VIP\freq_tuning_data_act.mat')   % activated group
data = inh;
NumCells = length(data);

%% multiplicative?

Multi = nan(NumCells,2);
for j = 1:NumCells
    x = data{j}.rate_tone;
    y = data{j}.rate_laser;
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    M0 = 1; %init c; c=1 no change
    
    [M,fval] = fminsearch(@(M)sum((M*x-y).^2),M0);     %multiplicative
    
    e = sqrt(fval); %return error
    Multi(j,1) = M;
    Multi(j,2) = e;
end

%% multiplicative?

Multi2 = nan(NumCells,2);
for j = 1:NumCells
    x = data{j}.rate_tone;
    y = data{j}.rate_laser;
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    M0 = 1; %init c; c=1 no change
    
    [M,fval] = fminsearch(@(M)sum((max(M*x,0)-y).^2),M0);     %multiplicative
    
    e = sqrt(fval); %return error
    Multi2(j,1) = M;
    Multi2(j,2) = e;
end

%% additive?

Add = nan(NumCells,2);
for j = 1:NumCells
    x = data{j}.rate_tone;
    y = data{j}.rate_laser;
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    M0 = 0; %init c; c=1 no change
    
    [A,fval2] = fminsearch(@(M)sum((M+x-y).^2),M0);       %additive
    e2 = sqrt(fval2); %return error
    Add(j,1) = A;
    Add(j,2) = e2;
end

%% additive?

Add2 = nan(NumCells,2);
for j = 1:NumCells
    x = data{j}.rate_tone;
    y = data{j}.rate_laser;
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    M0 = 0; %init c; c=1 no change
    
    [A,fval2] = fminsearch(@(M)sum((max(M+x,0)-y).^2),M0);       %additive
    e2 = sqrt(fval2); %return error
    Add2(j,1) = A;
    Add2(j,2) = e2;
end