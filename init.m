rng(42)
min = [0.024793 0.26 0.07 -0.7244 0.26 -0.08 0.4244 0 0 2 -15];
max = [0.024793 0.34 0.08 -0.4244 0.34 -0.07 0.7244 0 0 18 15];
par_path = [pwd '/data/para50.txt'];
log_path = [pwd '/data/log.txt'];
best_log = [pwd 'best_log.log'];
newgen = zeros(50,11);
for i = 1:50
    child = min + (max-min).*rand(1,11);
    newgen(i,:) = child;
end

p_format = '%f %f %f %f %f %f %f %f %f %f %f\n';
p_fid = fopen(par_path, 'w');
fprintf(p_fid,p_format,newgen.');
fclose(p_fid);

l_fid = fopen(log_path,'w');
fclose(l_fid);

for i = 1:50
    disp(i)
    %lastwarn('')
    
    par = newgen(i,:);
    [pts, self_cross] = evenpar(par);
    polyin = polyshape({pts(:,1 )},{pts(:,2)});
    [comx,comy] = centroid(polyin);
    s = num2str(i,' %03d');
    foil_path = [pwd '/data/ib' s '.inp'];
    [ib_fid msg] = fopen(foil_path, 'w');
    n = size(pts,1);
    % output to ib[i].inp
    % format: 1st line: npt; rest: x y coordinate
    
    if self_cross
        % write 0
        fprintf(ib_fid,'%d',0);
    else
        fprintf(ib_fid,'%d\n',n);
        fprintf(ib_fid,'%8.7f %8.7f\n',pts.');
    end
    fclose(ib_fid); 
end