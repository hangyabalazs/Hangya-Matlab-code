%% helper cell to plot info bar graph

load('C:\Balazs\_analysis\SOM\info5\allvar_ns.mat')
kl_som2
kl_nssom2=kl_som2
kl_nssom_err2=kl_som_err2
load('C:\Balazs\_analysis\SOM\info5\allvar_ws.mat')

%% Plot

figure;
bar([kl_pv kl_nssom2 kl_som kl_nso kl_wso])
hold on
errorbar([kl_pv kl_nssom2 kl_som kl_nso kl_wso],[kl_pv_err kl_nssom_err2 kl_som_err kl_ns_err kl_ws_err],'k+')
set(gca,'XTickLabel',{'PV HZout' 'NS SOM HZin' 'WS SOM HZout' 'NS HZout' 'WS HZout'})

%% Plot log2

figure;
bar(log2(exp(1))*[kl_pv kl_nssom2 kl_som kl_nso kl_wso])
hold on
errorbar(log2(exp(1))*[kl_pv kl_nssom2 kl_som kl_nso kl_wso],...
    log2(exp(1))*[kl_pv_err kl_nssom_err2 kl_som_err kl_ns_err kl_ws_err],'k+')
set(gca,'XTickLabel',{'PV HZout' 'NS SOM HZin' 'WS SOM HZout' 'NS HZout' 'WS HZout'})

%%

figure;     % HomeZoneOut
bar(log2(exp(1))*[kl_pv kl_som kl_nso kl_wso])
hold on
errorbar(log2(exp(1))*[kl_pv kl_som kl_nso kl_wso],...
    log2(exp(1))*[kl_pv_err kl_som_err kl_ns_err kl_ws_err],'k+')
title('HomeZoneOut')

figure;     % HomeZoneIn
bar(log2(exp(1))*[kl_pv2 kl_som2 kl_nso2 kl_wso2])
hold on
errorbar(log2(exp(1))*[kl_pv2 kl_som2 kl_nso2 kl_wso2],...
    log2(exp(1))*[kl_pv_err2 kl_som_err2 kl_ns_err2 kl_ws_err2],'k+')
title('HomeZoneIn')