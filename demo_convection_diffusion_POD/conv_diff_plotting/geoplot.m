function pmu = geoplot(probname,muval)
% -----rbMIT Software Copyright MIT 2006-07----
%  DBP  Huynh,NC Nguyen, AT Patera, G Rozza
% Henceforth, Software shall refer to the rbMIT_System software package: 
% the contents of the rbMIT_Library and rbMIT_Aux folders.
% 
% This License governs use of all accompanying Software, and your download and/or use of the Software constitutes acceptance 
% of this License.
% 
% You may use this Software for any non-commercial purpose, subject to the restrictions in this license. 
% Some purposes which can be non-commercial are teaching and academic research.
% 
% You may not use or distribute this Software or any derivative works in any form for commercial purposes. 
% Examples of commercial purposes would be licensing, leasing, or selling the Software, or distributing the 
% Software for use with commercial products.
% 
% You may modify this Software and distribute the modified Software for non-commercial purposes,
%  however, you may not grant rights to the Software or derivative works that are broader than those 
% provided by this License. For example, you may not distribute modifications of the Software under terms 
% that would permit commercial use.
% 
% In return, we require that you agree:
% 
% not to remove any copyright or other notices from the Software;
% 
% that if you distribute the Software in source or object form, you will include a verbatim copy of this license;
% 
% that if you distribute derivative works of the Software in source code form you do so only under a license
% that includes all of the provisions of this License, and if you distribute derivative works of the Software solely 
% in object form you do so only under a license that complies with this License;
% 
% that if you have modified the Software or created derivative works, and distribute such modifications or derivative
% works, you will cause the modified files to carry prominent notices so that recipients know that they are not receiving the original Software;
% 
% THAT THE SOFTWARE COMES "AS IS", WITH NO WARRANTIES: THIS MEANS NO EXPRESS, IMPLIED OR STATUTORY WARRANTY, 
% INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY WARRANTY 
% OF TITLE OR NON-INFRINGEMENT; YOU MUST PASS THIS DISCLAIMER ON WHENEVER YOU DISTRIBUTE THE SOFTWARE OR DERIVATIVE WORKS;
% 
% THAT NEITHER MIT NOR THE AUTHORS WILL BE LIABLE FOR ANY DAMAGES RELATED TO THE SOFTWARE OR THIS LICENSE, 
% INCLUDING DIRECT, INDIRECT, SPECIAL, CONSEQUENTIAL OR INCIDENTAL DAMAGES; YOU MUST PASS THIS LIMITATION OF LIABILITY 
% ON WHENEVER YOU DISTRIBUTE THE SOFTWARE OR DERIVATIVE WORKS;
% 
% that your rights under the License end automatically if you breach the License in any way;
% 
% that MIT reserves all rights not expressly granted to you in this License.
% 
% (Note: The wording of this license agreement is derived from a Microsoft Shared Source license agreement.)
% 
% License at http://augustine.mit.edu/rbMIT_SystemLicense.htm.
%
% rbMIT Software Part_II_and_IV
%
% ---------- 
% actually we only need Gaffs/Caffs and nodemap
load(strcat(probname,'_','SYMB'));

% geoplot
if (nargout == 0)
    copyfile(strcat(probname,'_','pdedom_save.m'),strcat(probname,'_','pdedom_save_tmp.m'))
    fid = fopen(strcat(probname,'_','pdedom_save_tmp.m'),'r+');
    for i = 1:12
        fgetl(fid);
    end
    fseek(fid,ftell(fid),-1);
    for i = 1:length(muval)
        fprintf(fid,'%s = %f;\n',strcat('mu',num2str(i)),muval(i));
    end
    fclose(fid);
    fclose('all');
    pdegplot(strcat(probname,'_','pdedom_save_tmp'));
    tmp = psym;
    for i = 1:P
        tmp = subs(tmp,strcat('mu',num2str(i)),muval(i));
    end
    prealtmp = tmp;
    hold on;
    for i = 1:numdoms
        tmp = sum(prealtmp(nodemap(i,:),:))/3;
        text(tmp(1),tmp(2),num2str(i));
    end
    axis equal; axis tight; hold off;
    delete(strcat(probname,'_','pdedom_save_tmp.m'));
end

% get new p coordinates
if (nargout > 0)
    % calculate new nodes locations and plot
    for i = 1:numdoms
        temp = Gaffs(i,:);
        for j = 1:P
            temp = subs(temp,strcat('mu',num2str(j)),muval(j));
        end
        Gaffsreal(i,:) = temp;
        temp = Caffs(i,:);
        for j = 1:P
            temp = subs(temp,strcat('mu',num2str(j)),muval(j));
        end
        Caffsreal(i,:) = temp;
    end
    ps = zeros(size(p,2),1);
    for i = 1:size(t,2)
        ps(t([1,2,3],i)) = t(4,i)*ones(3,1);
    end
    for i = 1:size(p,2)
        pmu(:,i) = [Gaffsreal(ps(i),1),Gaffsreal(ps(i),2);Gaffsreal(ps(i),3),Gaffsreal(ps(i),4)]*p(:,i)+[Caffsreal(ps(i),1);Caffsreal(ps(i),2)];
    end
end