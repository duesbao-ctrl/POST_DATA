function enabled = pd_catalog_enabled_mask(catalog, tableData)
%PD_CATALOG_ENABLED_MASK Evaluate parameter dependency conditions.

    enabled = true(numel(catalog), 1);
    for i = 1:numel(catalog)
        dependency = catalog(i).DependsOn;
        if isempty(dependency), continue; end
        parent = find(strcmp(dependency, {catalog.Name}), 1, 'first');
        if isempty(parent)
            error('postdata:MissingOptionDependency', ...
                'Option %s depends on missing option %s.', catalog(i).Name, dependency);
        end
        try
            parentValue = pd_parse_option_value(tableData{parent, 2}, catalog(parent).Type);
            enabled(i) = valuesEqual(parentValue, catalog(i).DependsValue);
        catch
            enabled(i) = false;
        end
    end
end

function equal = valuesEqual(actual, expected)
    if ischar(actual) || ischar(expected)
        equal = strcmpi(pd_to_char(actual), pd_to_char(expected));
    else
        equal = isequal(actual, expected);
    end
end
