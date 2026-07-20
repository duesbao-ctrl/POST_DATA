function pd_log(filePath, level, module, message)
%PD_LOG Append one diagnostic record. Logging failure never masks analysis.

    try
        folder = fileparts(filePath);
        if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
        rotateIfNeeded(filePath, 5 * 1024 * 1024);
        fid = fopen(filePath, 'a');
        if fid < 0, return; end
        cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
        timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS.FFF'); %#ok<DATST,TNOW1>
        cleanMessage = regexprep(pd_to_char(message), '[\r\n]+', ' ');
        fprintf(fid, '%s | %-7s | %-12s | %s\n', timestamp, upper(level), module, cleanMessage);
    catch
        % Diagnostics must not replace the original application error.
    end
end

function rotateIfNeeded(filePath, maximumBytes)
    info = dir(filePath);
    if isempty(info) || info(1).bytes <= maximumBytes, return; end
    archivePath = [filePath, '.1'];
    if exist(archivePath, 'file'), delete(archivePath); end
    movefile(filePath, archivePath, 'f');
end
