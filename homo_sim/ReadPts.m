function P = ReadPts(FileName)
    fid = fopen(FileName,'r');
    Npts = fscanf(fid,'%d',1);
    fscanf(fid,'%*[^\n]\n',1); % skip a line
    P = zeros(Npts,2); t = 1;
    while ~feof(fid)
        n = fscanf(fid,'%d%*[^\n]\n'); % read number of points
        for i = 1:n
            P(t,:) = fscanf(fid,'%f%f',[1,2]); % read 2 numbers
            fscanf(fid,'%*[^\n]\n',1); % skip a line
            t = t + 1;
        end
    end
    fclose(fid);
end