function thrombin
read_file
end
function read_file
fid = fopen('equations.txt');
fout = fopen('edge_list.txt','w');
tab=';';
tline1 = fgetl(fid);
tline2 = fgetl(fid);
j=0;vertex_list={};
while ischar(tline1)    
    left=tline1(2:end);
    ind=find(ismember(vertex_list,left))
    if isempty(ind)
        j=j+1;
    vertex_list{j}=left;
    end
    openid=strfind(tline2,'[');
    closingid=strfind(tline2,']');
    for i=1:length(openid)
        right=tline2(openid(i):closingid(i));
        
    ind=find(ismember(vertex_list,right));
    
    if isempty(ind) 
        j=j+1;
    vertex_list{j}=right;
    end
        fprintf(fout,['%s',tab,'%s\n'],left,right);
        idleft=find(ismember(vertex_list,left))
        idright=find(ismember(vertex_list,right))
        A(idleft,idright)=1;
    end
    
tline1 = fgetl(fid);
tline2 = fgetl(fid);
end
for ii=1:34
    tree=bfs(A,ii,vertex_list)
    sizetree(ii)=length(tree);
end
fclose(fid);
fclose(fout);
end
function tree=bfs(A,s,vertex_list)
%mark,tree
n=max(size(A))
visited=zeros(n,1);
mark(1)=s;
tree(1)=s;
visited(s)=1;
while length(mark)>0
    i=mark(1);
    
    tem=find(A(i,:)>0);
    if length(tem)>0
        for j=1:length(tem)
            if visited(tem(j))==0
                mark(end+1)=tem(j);
                tree(end+1)=tem(j);
                visited(tem(j))=1;
            end
        end
    end
    mark(1)=[];
end
end