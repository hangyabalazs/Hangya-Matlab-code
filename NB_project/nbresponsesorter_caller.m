loadcb
I = CELLIDLIST(3444:end);
selstr = ['"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
    'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
    'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
I2 = selectcell(selstr);
I3 = I2(1143:end);
error_list = nbresponsesorter(I3,true);