% from Hyun

%% load fitted version

% load('C:\Balazs\_analysis\VIP\revision_figures\fit_variables_inh.mat')
load('C:\Balazs\_analysis\VIP\revision_figures\fit_variables_act.mat')
NumCells = length(fitted_nostim);
data1 = fitted_nostim;
data2 = fitted_stim;

%% multiplicative?

Multi = nan(NumCells,2);
for j = 1:NumCells
    x = data1{j};
    y = data2{j};
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    M0 = 1; %init c; c=1 no change
    
    [M,fval] = fminsearch(@(M)sum((M*x-y).^2),M0);     %multiplicative
    
    e = sqrt(fval); %return error
    Multi(j,1) = M;
    Multi(j,2) = e;
end

%% additive?

Add = nan(NumCells,2);
for j = 1:NumCells
    x = data1{j};
    y = data2{j};
    x = x(~isnan(x));
    y = y(~isnan(y));
    
    M0 = 0; %init c; c=1 no change
    
    [A,fval2] = fminsearch(@(M)sum((M+x-y).^2),M0);       %additive
    e2 = sqrt(fval2); %return error
    Add(j,1) = A;
    Add(j,2) = e2;
end