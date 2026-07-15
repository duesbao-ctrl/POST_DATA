function value = pd_ui_text(key, varargin)
%PD_UI_TEXT Return a localized UI string without non-ASCII source literals.
% MATLAB R2016b on Windows reads M-files with the active system code page.
% UI translations therefore live in an explicitly decoded UTF-8 resource.
% The English key is also the safe fallback when the resource is unavailable.

    key = pd_to_char(key);
    persistent catalog
    if isempty(catalog)
        catalog = loadCatalog();
    end

    if isKey(catalog, key)
        value = catalog(key);
    else
        value = key;
    end
    if ~isempty(varargin)
        value = sprintf(value, varargin{:});
    end
end

function catalog = loadCatalog()
    catalog = containers.Map('KeyType', 'char', 'ValueType', 'char');
    resourcePath = fullfile(fileparts(mfilename('fullpath')), ...
        'resources', 'ui_zh_CN.tsv');
    fid = fopen(resourcePath, 'r', 'n', 'UTF-8');
    if fid < 0
        warning('postdata:uiText:MissingResource', ...
            'UI translation resource is unavailable: %s', resourcePath);
        return;
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lineNumber = 0;
    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        lineNumber = lineNumber + 1;
        if lineNumber == 1 && ~isempty(line) && double(line(1)) == 65279
            line = line(2:end);
        end
        if isempty(line) || line(1) == '#'
            continue;
        end
        separator = find(line == char(9), 1, 'first');
        if isempty(separator)
            error('postdata:uiText:MalformedResource', ...
                'Malformed UI translation at line %d in %s.', ...
                lineNumber, resourcePath);
        end
        source = strtrim(line(1:separator - 1));
        translation = line(separator + 1:end);
        if isempty(source)
            error('postdata:uiText:EmptyKey', ...
                'Empty UI translation key at line %d in %s.', ...
                lineNumber, resourcePath);
        end
        if isKey(catalog, source)
            error('postdata:uiText:DuplicateKey', ...
                'Duplicate UI translation key "%s" at line %d.', ...
                source, lineNumber);
        end
        catalog(source) = translation;
    end
end
