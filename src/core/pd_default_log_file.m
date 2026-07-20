function filePath = pd_default_log_file()
%PD_DEFAULT_LOG_FILE Return the per-user POST_DATA2 log path.

    logDir = fullfile(prefdir, 'POST_DATA2', 'logs');
    if ~exist(logDir, 'dir')
        mkdir(logDir);
    end
    filePath = fullfile(logDir, 'postdata.log');
end
