function manifest = pd_request_manifest(request)
%PD_REQUEST_MANIFEST Remove runtime callbacks from reproducibility metadata.

    manifest = request;
    if isfield(manifest, 'execution')
        manifest.execution.cancelCallback = '<runtime callback>';
        manifest.execution.progressCallback = '<runtime callback>';
    end
end
