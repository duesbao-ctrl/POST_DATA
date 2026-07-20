function pd_execution_checkpoint(request, fraction, message)
%PD_EXECUTION_CHECKPOINT Report progress and honor cooperative cancellation.

    drawnow;
    if request.execution.cancelCallback()
        error('postdata:UserCancelled', 'Analysis cancelled by user.');
    end
    request.execution.progressCallback(max(0, min(1, fraction)), message);
end
