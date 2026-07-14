function pd_write_numeric_csv(filePath, headers, data)
%PD_WRITE_NUMERIC_CSV Write a numeric matrix with a CSV header.

    fid = fopen(filePath, 'w');
    if fid < 0, error('postdata:CsvOpenFailed', 'Cannot create CSV: %s', filePath); end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n', strjoin(headers, ','));
    if isempty(data), return; end
    format = [repmat('%.16g,', 1, size(data, 2) - 1), '%.16g\n'];
    fprintf(fid, format, data.');
end
