function phi = loadUGmex(file, plot)

    [phi, bb] = loadug(file);

    if nargin > 1 && plot == 1
        [XX,YY,ZZ] = meshgrid(bb(:,1), bb(:,2), bb(:,3));
        slice(XX, YY, ZZ, phi, (bb(1,1) + bb(end,1)) / 2, ...
                               (bb(1,2) + bb(end,2)) / 2, ...
                               (bb(1,3) + bb(end,3)) / 2);
        xlabel('X'); ylabel('Y'); zlabel('Z');
                   
        figure;
        isosurface(XX, YY, ZZ, phi, 0);
    end
    
end