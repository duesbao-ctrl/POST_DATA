function fullPath = get_slurm_txt_fullpath(dirPath)
%GET_SLURM_TXT_FULLPATH Get full path of unique slurm* file in a folder.
%   fullPath = GET_SLURM_TXT_FULLPATH(dirPath)
%
% Input:
%   dirPath : target directory.
%
% Output:
%   fullPath: full file path (char).
%
% Rules:
%   - exactly one match for pattern slurm* is required
%   - if zero or multiple matches, throw error

    if nargin < 1 || isempty(dirPath)
        error('get_slurm_txt_fullpath:MissingInput', 'dirPath is required.');
    end
    if isstring(dirPath)
        dirPath = char(dirPath);
    end
    if ~ischar(dirPath)
        error('get_slurm_txt_fullpath:InvalidInputType', 'dirPath must be char or string.');
    end
    if ~exist(dirPath, 'dir')
        error('get_slurm_txt_fullpath:DirNotFound', 'Directory not found: %s', dirPath);
    end

    files = dir(fullfile(dirPath, 'slurm*'));
    files = files(~[files.isdir]);

    if isempty(files)
        error('get_slurm_txt_fullpath:NoMatch', ...
            'No file matches pattern slurm* in: %s', dirPath);
    end
    if numel(files) > 1
        names = {files.name};
        error('get_slurm_txt_fullpath:MultipleMatches', ...
            'Multiple slurm* files found in %s: %s', dirPath, strjoin(names, ', '));
    end

    fullPath = fullfile(dirPath, files(1).name);
end
