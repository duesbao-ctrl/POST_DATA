function rootDir = postdata_startup()
%POSTDATA_STARTUP Add POST_DATA source folders to the MATLAB path.
% MATLAB R2016b compatible. No package folders are used.

    rootDir = fileparts(mfilename('fullpath'));
    sourceDir = fullfile(rootDir, 'src');
    if ~exist(sourceDir, 'dir')
        error('postdata:MissingSourceDirectory', ...
            'Source directory not found: %s', sourceDir);
    end
    addpath(genpath(sourceDir));
end
