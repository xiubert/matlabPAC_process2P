function [d,c] = testSandbox(name)
%testSandbox  Fresh scratch folder for a test, deleted when c goes out of scope.
%
%   [d,c] = testSandbox(name)
%
%   Returns an empty folder under tempdir named for the test, plus an
%   onCleanup that removes it. Hold onto c:
%
%       [tmp,cleanup] = testSandbox('myTest');   %#ok<ASGLU>
%
%   Tests write synthesized animal folders and redirected save targets here so
%   they never touch the recorded fixtures. Any pre-existing folder of the
%   same name is cleared first, so a test that died mid-run last time does not
%   leave state behind for the next one.
%
%   See also testConfig, requireFixture.

d = fullfile(tempdir,['matlabPACtest_' char(name)]);
if isfolder(d); rmdir(d,'s'); end
mkdir(d);
c = onCleanup(@() cleanupDir(d));
end

function cleanupDir(d)
if isfolder(d)
    try, rmdir(d,'s'); catch, end   % a locked .mat must not fail the test
end
end
