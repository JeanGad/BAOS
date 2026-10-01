function write_csv(fname, header, M)
% WRITE_CSV  Toolbox-free CSV writer (numeric matrix + header cellstr).
fid = fopen(fname, 'w'); fprintf(fid, '%s', header{1});
fprintf(fid, ',%s', header{2:end}); fprintf(fid, '\n');
fmt = [repmat('%.10g,', 1, size(M,2)-1) '%.10g\n'];
fprintf(fid, fmt, M.'); fclose(fid);
end
