function myfunc(gen)
% disp(gen)
log_fid = fopen('test.txt','a');

best_par = [1 2 3 4 5 6 7 8 9 10 11];
fprintf(log_fid,'%d %f %f %f %f %f %f %f %f %f %f %f\n',[gen best_par]);
sprintf('%d %f %f %f %f %f %f %f %f %f %f %f\n',[gen best_par])
