function info = pd_version()
%VERSION Return POST_DATA software and result-schema version information.
% MATLAB R2016b compatible.

    info = struct();
    info.name = 'POST_DATA';
    info.version = '0.6.0';
    info.resultSchemaVersion = 1;
    info.minimumMatlabRelease = 'R2016b';
end
