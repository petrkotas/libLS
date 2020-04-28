function [grid, data] = readArray(file, plot)

    fd = fopen(file, 'rb');
    
    grid.dim  = fread(fd, 1, 'int');
    grid.dims = fread(fd, 3, 'int');
    grid.dx   = fread(fd, 3, 'double');
    grid.bl   = fread(fd, 3, 'double');
    grid.tr   = fread(fd, 3, 'double');
    
    num_gridpts = grid.dims(1) * grid.dims(2) * grid.dims(3);    
    data = fread(fd, num_gridpts, '*double');

    fclose(fd);
   
    data = reshape(data, grid.dims');
    data = permute(data, [2,1,3]);
    
    if nargin > 1 && plot == 1
        [XX,YY,ZZ] = meshgrid(grid.bl(1):grid.dx(1):grid.tr(1), ...
                              grid.bl(2):grid.dx(2):grid.tr(2), ...
                              grid.bl(3):grid.dx(3):grid.tr(3));
        figure;
        xslice = [(2*grid.bl(1) + grid.tr(1))/3, (grid.bl(1) + 2*grid.tr(1))/3];
        yslice = [(2*grid.bl(2) + grid.tr(2))/3, (grid.bl(2) + 2*grid.tr(2))/3];
        zslice = [(2*grid.bl(3) + grid.tr(3))/3, (grid.bl(3) + 2*grid.tr(3))/3];
        h = slice(XX, YY, ZZ, data, xslice, ...
                                    yslice, ...
                                    zslice);

        
        xlabel('X'); ylabel('Y'); zlabel('Z');
                   
        figure;
        isosurface(data, 0);
    end
    
end