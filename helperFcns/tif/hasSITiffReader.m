function [tf,apiVer] = hasSITiffReader(pkgParent)
%HASSITIFFREADER  True when Vidrio's ScanImageTiffReader package is usable.
%
%   tf            = hasSITiffReader
%   [tf,apiVer]   = hasSITiffReader(pkgParent)
%
%   PKGPARENT is the folder that *contains* +ScanImageTiffReader; the package
%   folder itself is never added to the path (MATLAB resolves packages from
%   the parent).  Defaults to /home/pac/Documents/MATLAB.
%
%   The mex is actually exercised through apiVersion(), so a
%   present-but-unloadable mex (wrong platform, missing libScanImageTiffReader
%   shared library, no LD_LIBRARY_PATH) reports false here instead of
%   throwing halfway through a read.
%
%   See also flattenTif, readSCIMtif

if nargin < 1 || isempty(pkgParent)
    pkgParent = '/home/pac/Documents/MATLAB';
end
pkgParent = char(pkgParent);

%add the parent folder if the package lives there and it isn't on the path yet
if isfolder(fullfile(pkgParent,'+ScanImageTiffReader')) && ...
        ~any(strcmp(pkgParent,strsplit(path,pathsep)))
    addpath(pkgParent);
end

try
    apiVer = ScanImageTiffReader.ScanImageTiffReader.apiVersion();
    tf     = true;
catch
    apiVer = '';
    tf     = false;
end

end