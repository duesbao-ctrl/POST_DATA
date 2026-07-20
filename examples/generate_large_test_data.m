function paths = generate_large_test_data(outputDir)
%GENERATE_LARGE_TEST_DATA Create deterministic 1D, 2D, cluster, and mass-v fixtures.
% Compatible with MATLAB R2016b. Files contain three timesteps each.

    if nargin < 1 || isempty(outputDir)
        rootDir = fileparts(fileparts(mfilename('fullpath')));
        outputDir = fullfile(rootDir, 'fixtures', 'generated');
    end
    if ~exist(outputDir, 'dir')
        [ok, message] = mkdir(outputDir);
        if ~ok, error('postdata:FixtureDirectoryFailed', '%s', message); end
    end

    paths = struct();
    paths.chunk1d = fullfile(outputDir, 'large_bin1d_dx_0.025.txt');
    paths.chunk2d = fullfile(outputDir, ...
        'large_bin2d_dx_0.05_dy_0.06_Lz_1.txt');
    paths.cluster = fullfile(outputDir, 'large_cluster_chunk.txt');
    paths.massv = fullfile(outputDir, 'large_mass_v.txt');
    writeChunk1d(paths.chunk1d);
    writeChunk2d(paths.chunk2d);
    writeCluster(paths.cluster);
    writeMassV(paths.massv);
end

function writeMassV(filePath)
    fid = openOutput(filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '# generated large mass-v fixture\n');
    fprintf(fid, '# 501 velocity bins and 3 timesteps\n');
    fprintf(fid, ['# Chunk mass1ArealDensity mass2ArealDensity ', ...
        'massArealDensity\n']);
    velocity = (-5000:20:5000).';
    timesteps = [100, 200, 300];
    for s = 1:numel(timesteps)
        center1 = 400 + 350 * (s - 1);
        center2 = -900 + 220 * (s - 1);
        mass1 = 0.002 + 0.030 .* exp(-((velocity - center1) ./ 1050).^2);
        mass2 = 0.001 + 0.020 .* exp(-((velocity - center2) ./ 1450).^2);
        fprintf(fid, '%d %d %d\n', timesteps(s), numel(velocity), 0);
        for i = 1:numel(velocity)
            fprintf(fid, '%.8g %.8g %.8g %.8g\n', velocity(i), ...
                mass1(i), mass2(i), mass1(i) + mass2(i));
        end
    end
end

function writeChunk1d(filePath)
    fid = openOutput(filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '# generated large 1D fixture\n');
    fprintf(fid, '# 401 bins and 3 timesteps\n');
    fprintf(fid, '# Chunk Coord1 Ncount c_rho\n');
    n = 401;
    x = (0:n - 1).' .* 0.025;
    timesteps = [100, 200, 300];
    for s = 1:numel(timesteps)
        center = 3.2 + 0.8 * (s - 1);
        rho = 0.8 + 1.7 .* exp(-((x - center) ./ 1.25).^2) + ...
            0.15 .* sin(2 .* pi .* x ./ 2.5 + 0.4 .* s);
        ncount = max(0, round(8 + 32 .* exp(-((x - center) ./ 1.35).^2) + ...
            3 .* sin(2 .* pi .* x ./ 2.5 + 0.4 .* s)));
        fprintf(fid, '%d %d %d\n', timesteps(s), n, sum(ncount));
        for i = 1:n
            fprintf(fid, '%d %.8g %d %.8g\n', ...
                i, x(i), ncount(i), rho(i));
        end
    end
end

function writeChunk2d(filePath)
    fid = openOutput(filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '# generated large 2D fixture\n');
    fprintf(fid, '# 61 x 41 grid and 3 timesteps\n');
    fprintf(fid, '# Chunk Coord1 Coord2 Ncount c_rho\n');
    nx = 61;
    ny = 41;
    xValues = (0:nx - 1) .* 0.05;
    yValues = (0:ny - 1) .* 0.06;
    n = nx * ny;
    timesteps = [100, 200, 300];
    for s = 1:numel(timesteps)
        fprintf(fid, '%d %d %d\n', timesteps(s), n, 5 * n);
        row = 0;
        for iy = 1:ny
            y = yValues(iy);
            for ix = 1:nx
                x = xValues(ix);
                row = row + 1;
                movingX = 0.85 + 0.18 * (s - 1);
                pore1 = ((x - movingX) / 0.48)^2 + ((y - 1.10) / 0.55)^2 < 1;
                pore2 = ((x - 2.15) / 0.38)^2 + ((y - 1.72) / 0.42)^2 < 1;
                channel = abs(y - (0.45 + 0.20 * sin(2.4 * x + 0.3 * s))) < 0.075;
                if pore1 || pore2 || channel
                    ncount = 0;
                else
                    ncount = 2 + mod(ix + 2 * iy + s, 4);
                end
                rho = 0.25 + 0.32 * ncount + 0.12 * sin(1.7 * x) * cos(1.4 * y);
                fprintf(fid, '%d %.8g %.8g %d %.8g\n', ...
                    row, x, y, ncount, rho);
            end
        end
    end
end

function writeCluster(filePath)
    fid = openOutput(filePath);
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '# generated large cluster fixture\n');
    fprintf(fid, '# 1200 clusters and 3 timesteps\n');
    fprintf(fid, '# Chunk c_x c_y c_z vx vy vz Ncount\n');
    n = 1200;
    timesteps = [100, 200, 300];
    index = (1:n).';
    for s = 1:numel(timesteps)
        x = 12 .* mod(index .* 0.61803398875 + 0.07 * s, 1);
        y = 8 .* mod(index .* 0.41421356237 + 0.11 * s, 1);
        z = 3 .* mod(index .* 0.73205080757 + 0.05 * s, 1);
        vx = 1.8 .* sin(index .* 0.031 + 0.2 * s);
        vy = 1.1 .* cos(index .* 0.027 + 0.3 * s);
        vz = 0.5 .* sin(index .* 0.019);
        ncount = 2 + mod(index .* 17 + index .* index + 13 * s, 500);
        fprintf(fid, '%d %d %d\n', timesteps(s), n, sum(ncount));
        for i = 1:n
            fprintf(fid, '%d %.8g %.8g %.8g %.8g %.8g %.8g %d\n', ...
                i, x(i), y(i), z(i), vx(i), vy(i), vz(i), ncount(i));
        end
    end
end

function fid = openOutput(filePath)
    fid = fopen(filePath, 'w');
    if fid < 0
        error('postdata:FixtureOpenFailed', 'Cannot create fixture: %s', filePath);
    end
end
