function p = requireFixture(p,label)
%requireFixture  Skip the calling test when a recorded fixture is absent.
%
%   p = requireFixture(p,label)
%
%   Throws matlabPAC:test:fixtureMissing if p is neither a file nor a folder.
%   runTests treats that identifier as SKIP rather than FAIL, so the suite
%   still reports a clean result on a machine that has the repo but not the
%   lab data drive. Any other error from a test is a real failure.
%
%   label is what shows up in the skip message, e.g. 'TO0003 animal folder'.
%
%   See also testConfig, runTests.

if nargin < 2; label = 'fixture'; end
p = char(p);
if ~isfolder(p) && ~isfile(p)
    error('matlabPAC:test:fixtureMissing','%s not found: %s', label, p);
end
end
