function app = postdata_app()
%POSTDATA_APP Launch the MATLAB R2016b-compatible POST_DATA application.
%
%   postdata_app
%   app = postdata_app()

    postdata_startup();
    existing = findall(0, 'Type', 'figure', 'Tag', 'POST_DATA2_MainFigure');
    if isempty(existing)
        existing = findall(0, 'Type', 'figure', 'Name', 'POST_DATA2 MATLAB 后处理软件');
    end
    if ~isempty(existing)
        delete(existing);
        drawnow;
    end
    app = PostDataApp();
    setappdata(app.Figure, 'POST_DATA_App', app);
    if nargout == 0
        clear app;
    end
end
