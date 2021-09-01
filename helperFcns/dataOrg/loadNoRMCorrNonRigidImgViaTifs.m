function moCorrImgNonRigid = loadNoRMCorrNonRigidImgViaTifs(dataPath,varargin)
if ~isempty(varargin)
    Ftifs = varargin{1};
    [Ftifs.folder] = deal([dataPath filesep 'NoRMCorred']);
    normCorreTifs = strrep({Ftifs.name},'.tif','_NoRMCorre.tif')';
    [Ftifs.name] = normCorreTifs{:};
    NoRMCorreTifs = Ftifs;
else
    NoRMCorreTifs = dir([dataPath filesep 'NoRMCorred' filesep '*_NoRMCorre.tif']);
end
[Ycon,~] = concatenate_files(NoRMCorreTifs);
moCorrImgNonRigid = single(Ycon); % convert to single precision
moCorrImgNonRigid = moCorrImgNonRigid - min(moCorrImgNonRigid(:));
end