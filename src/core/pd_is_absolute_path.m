function tf = pd_is_absolute_path(pathValue)
%PD_IS_ABSOLUTE_PATH True for Windows drive/UNC and Unix absolute paths.

    pathValue = pd_to_char(pathValue);
    tf = ~isempty(regexp(pathValue, ...
        '^(?:[A-Za-z]:[\\/]|[\\/]{2}|/)', 'once'));
end
