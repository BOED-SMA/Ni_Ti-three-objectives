function [PF, PF_indices] = Find_pareto_front_multi2(objs)

% author: Shahin Boluki
% date: Oct 2017
% modified to add an additional output PF

%Find pareto front
%objs : n*m n=number of samples m=number of objectives
%works for any nymber of objectives
%Assumption: minimizing
PF_indices = [];
PF = [];
n = size(objs,1);
m=size(objs,2);
for i=1:1:n
    better_set=cell(1,m);
    better_equal_set=cell(1,m);
    for j=1:1:m
        better_set{j}=find(objs(:,j)<objs(i,j));
        better_equal_set{j}=find(objs(:,j)<=objs(i,j));
    end
    flag_not_pareto=0;
    for j=1:1:m
        
        temp_inds = better_set{j};
        for kk=1:1:m
            if kk==j
                continue;
            end
            if isempty(temp_inds)
               break; 
            end
            temp_inds=intersect(temp_inds,better_equal_set{kk});
        end
        if(~isempty(temp_inds))
            flag_not_pareto = 1;
            break;
        end

    end
    if flag_not_pareto==0
        PF_indices = [PF_indices;i];
        PF = [PF; objs(i, :)];
    end
end

end

