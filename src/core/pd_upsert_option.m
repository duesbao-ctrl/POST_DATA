function options = pd_upsert_option(options, name, value)
%PD_UPSERT_OPTION Replace or append one name-value option.

    keep = true(size(options));
    for i = 1:2:numel(options)
        key = options{i};
        if pd_is_text_scalar(key) && strcmpi(pd_to_char(key), name)
            keep(i:min(i + 1, numel(options))) = false;
        end
    end
    options = options(keep);
    options = [options, {name, value}];
end
