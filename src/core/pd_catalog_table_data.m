function data = pd_catalog_table_data(catalog)
%PD_CATALOG_TABLE_DATA Convert option metadata into GUI table data.

    data = cell(numel(catalog), 3);
    for i = 1:numel(catalog)
        data{i, 1} = catalog(i).Name;
        data{i, 2} = pd_format_option_value(catalog(i).Default);
        data{i, 3} = catalog(i).Description;
    end
end
