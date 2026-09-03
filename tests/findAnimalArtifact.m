function f = findAnimalArtifact(dataPath,name)
%findAnimalArtifact  Locate an animal artifact under either layout.
%
%   f = findAnimalArtifact(dataPath,name)
%
%   Artifacts used to sit loose in the animal folder and now live under
%   analysis/<run>/ (see animalPaths). A fixture path hardcoded to either one
%   breaks when a folder is migrated, so tests resolve through this instead:
%   the flat location first, then any run, newest run last so it wins.
%
%   Returns '' when the artifact is nowhere, which lets a test decide between
%   skipping and failing.
%
%   See also animalPaths, migrateAnimalArtifacts, requireFixture.

f = '';
flat = fullfile(dataPath,name);
if isfile(flat); f = flat; return, end

d = dir(fullfile(dataPath,'analysis','*',name));
if isempty(d); return, end
[~,ord] = sort([d.datenum]);
d = d(ord(end));
f = fullfile(d.folder,d.name);
end
