function value = pd_get_analysis_option(options, name, defaultValue)
%PD_GET_ANALYSIS_OPTION Read one case-insensitive name-value option.

    value = defaultValue;
    for i = 1:2:numel(options)
        if pd_is_text_scalar(options{i}) && strcmpi(pd_to_char(options{i}), name)
            value = options{i + 1};
            return;
        end
    end
end
