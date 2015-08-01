setpref('Internet','E_mail','hangyab@koki.hu');
% setpref('Internet', 'SMTP_Server', 'smtp.mail.yahoo.com');
setpref('Internet', 'SMTP_Server', 'smtp.koki.hu');
sendmail('balazs.cshl@gmail.com', 'Hello From MATLAB!');