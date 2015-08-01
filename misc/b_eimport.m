function [data_control,data_patient] = b_eimport
%EIMPORT   Import 'Ellsticity' data.
%   [D1,D2] = EIMPORT imports control (D1) and patient (D2) group data.

% Input argument check
error(nargchk(0,0,nargin))

% Import
where = 'F:\balazs\MersichBea2\';
fn_control = [where 'kontroll.xls'];
fn_patient = [where 'patient.xls'];
data_control = xlsread(fn_control);
data_patient = xlsread(fn_patient);