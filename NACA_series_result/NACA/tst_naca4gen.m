%first  = [2000 4000 6000];
%second = [200 400 600];
last_two = [12 20 30 50]; %thin,mid,thicc,too thicc
first_two = [0000 2200 2400 2600 4200 4400 4600 6200 6400 6600]
% last_two = 12;
% first_two = 0

%numberlist

logfile = fopen([pwd '/foil/naca_gen.log'],'w')
for num12 = first_two
    for num34 = last_two
        s = num2str(num12 + num34,'%04d');
        disp(s)
        filename = [pwd '/foil/NACA' s '.inp'];
        %         filename = [pwd 'degug.inp'];
        [fileID,message] = fopen(filename,'w');
        iaf.designation=s;
        % designation='0008';
        iaf.n=30;
        iaf.HalfCosineSpacing=1;
        iaf.wantFile=0;
        iaf.datFilePath='./'; % Current folder
        iaf.is_finiteTE=0;
        
        af = naca4gen(iaf);
        
        % plot(af.x,af.z,'bo-')
        
        % plot(af.xU,af.zU,'bo-')
        % hold on
        % plot(af.xL,af.zL,'ro-')
        
        t=0.015;
        pt = interparc(101, af.x, af.z);
        
        targetd = 0.015;
        d = pdist([pt(50,:);pt(51,:)]);
        while 1
            newn = round(100*d/targetd);
%             disp(newn);
            newpt = interparc(newn, af.x, af.z);
            newpt(newn,:) = []
            newpt(:,1) = newpt(:,1)-0.5
%             plot(newpt(:,1),newpt(:,2),'bo')
%             hold on
%             axis equal
            d = pdist([newpt(1,:);newpt(2,:)])
            if abs((d - targetd)/targetd) < 1
                break
            else
                fprintf(logfile, '  d = %f\n, regen', d)
            end
        end
        
        fprintf(fileID,'%d\n',newn-1)
        fprintf(fileID,'%f %f\n',newpt.');
        fclose(fileID);
        fprintf(logfile, 'NACA %s done!\n', s)
    end
end
fclose(logfile)