function phi = loadUGtxt(file, plot)

    fd = fopen(file, 'rt');
    C  = textscan(fd, '%s', 'Delimiter', '\r\n');
    fclose(fd);
    
    dx      = sscanf(C{1}{3}, '%f');
    maxX(1) = sscanf(C{1}{5}, '%*s %f');
    minX(1) = sscanf(C{1}{6}, '%*s %f');
    
    maxX(2) = sscanf(C{1}{7}, '%*s %f');
    minX(2) = sscanf(C{1}{8}, '%*s %f');
    
    maxX(3) = sscanf(C{1}{9}, '%*s %f');
    minX(3) = sscanf(C{1}{10}, '%*s %f');
    
    count = sscanf(C{1}{11}, '%*s %f');
    
    phi = cellfun(@(c)sscanf(c, '%f'), C{1}(12:end-1));
    
    disp(length(phi));
    disp(count);

    M = int32((maxX(1) - minX(1)) / dx);
    N = int32((maxX(2) - minX(2)) / dx);
    P = int32((maxX(3) - minX(3)) / dx);
    
    disp([M, N, P]);
    disp(M*N*P);

    X = linspace(minX(1), maxX(1), M);
    Y = linspace(minX(2), maxX(2), N);
    Z = linspace(minX(3), maxX(3), P);
    
    phi = reshape(phi, [length(Y), length(X), length(Z)]);

    
    if nargin > 1 && plot == 1
        [XX,YY,ZZ] = meshgrid(X, Y, Z);
        slice(XX, YY, ZZ, phi, (minX(1) + maxX(1)) / 2, ...
                               (minX(2) + maxX(2)) / 2, ...
                               (minX(3) + maxX(3)) / 2);
        xlabel('X'); ylabel('Y'); zlabel('Z');
                   
        figure;
        isosurface(XX, YY, ZZ, phi, 0);
    end
    
end