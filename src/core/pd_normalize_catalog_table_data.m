function data = pd_normalize_catalog_table_data(catalog, data)
%PD_NORMALIZE_CATALOG_TABLE_DATA Canonicalize editable option-cell display.

    if size(data, 1) ~= numel(catalog)
        error('postdata:OptionRowMismatch', ...
            'Option table does not match the active catalog.');
    end
    for i = 1:numel(catalog)
        data{i, 2} = pd_format_option_editor_value(data{i, 2}, ...
            catalog(i).Type);
    end
end
