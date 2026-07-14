function figures = pd_style_result_figures(result, options)
%PD_STYLE_RESULT_FIGURES Style every standalone figure contained in a result.

    if nargin < 2, options = struct(); end
    figures = collectFigures(result);
    for i = 1:numel(figures)
        axesFound = findall(figures(i), 'Type', 'axes');
        for j = 1:numel(axesFound)
            tag = get(axesFound(j), 'Tag');
            if any(strcmpi(tag, {'legend','Colorbar'})), continue; end
            pd_apply_publication_style(axesFound(j), options);
        end
        set(figures(i), 'Color', 'w');
    end
end

function figures = collectFigures(value)
    figures = gobjects(0, 1);
    if ~isstruct(value), return; end
    names = fieldnames(value);
    for k = 1:numel(value)
        for i = 1:numel(names)
            fieldValue = value(k).(names{i});
            lowerName = lower(names{i});
            if (~isempty(strfind(lowerName, 'fig')) || strcmp(lowerName, 'figure')) && ... %#ok<STREMP>
                    ~isempty(fieldValue)
                candidates = fieldValue(:);
                for j = 1:numel(candidates)
                    if ishghandle(candidates(j)) && ...
                            strcmp(get(candidates(j), 'Type'), 'figure')
                        figures(end + 1, 1) = candidates(j); %#ok<AGROW>
                    end
                end
            elseif isstruct(fieldValue)
                nested = collectFigures(fieldValue);
                for j = 1:numel(nested)
                    if ~any(figures == nested(j))
                        figures(end + 1, 1) = nested(j); %#ok<AGROW>
                    end
                end
            end
        end
    end
end
