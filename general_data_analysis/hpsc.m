%% Peri-stimulus count

fn = 'X:\In_Vivo\_analysis\acsady_csoport\all cells PSTH_counts.xls';
x = xlsread(fn,'sheet2','AW1:AW7000');
y = xlsread(fn,'sheet2','AX1:AX7000');

x = x * 1000;
y = y * 1000;
r1 = 0;
r2 = 0;
for k = 1:20
    lg1 = (x > y(k) - 20) & (x < y(k));
    lg2 = (x > y(k)) & (x < y(k) + 20);
    r1 = r1 + length(find(lg1));   % before
    r2 = r2 + length(find(lg2));   % after
end