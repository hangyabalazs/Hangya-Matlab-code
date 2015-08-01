function preallocate(xdim,ydim,v,varargin)
%PREALLOCATE   Preallocate memory space for variables.
%   PREALLOCATE(XD,YD,X,V1,V2,V3,...) preallocates memory space for
%   variables V1, V2, V3,... of size XD-by-YD with spacing values passed in
%   X (typically 0, 1, or NaN) in the workspace of the caller function.
%
%   See also ZEROS, ONES and NAN.