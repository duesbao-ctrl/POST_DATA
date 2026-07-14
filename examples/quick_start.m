% QUICK_START Minimal POST_DATA example for MATLAB R2016b and later.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(rootDir);
postdata_startup();

request = pd_create_request('network2d');
request.baseDir = fullfile(rootDir, 'fixtures', 'sample');
request.filePath = 'bin2d_dx_0.5_dy_0.5_Lz_1.txt';
request.analysisOptions = {'ThresholdN', 1, 'GeometryMode', 'original'};
request.makePlots = false;

result = postdata_run(request);
fprintf('Porosity: %.6f\n', result.global.porosity);

figure('Color', 'w');
pd_render_result(gca, result);
