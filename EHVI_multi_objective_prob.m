function eh_integral = EHVI_multi_objective_prob(mu, sigma, front_array, reference_point)
% A function to compute expected hyper-volumn improvement (EHVI) of m 
%objectives (m >= 2). The result is the EHVI of a single point with reference 
%to an approximation set. The objectives are assumed to be Gaussian
%distributed. The function follows the algorithm discribed in reference: 
%http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5949880
%
% eh_integral: the value of EHVI of one point
%
% mu, sigma: the mean and standard variance of the objectives which stored
%in a vector of length m. 
%
% front_array: an n*m matrix stored the data of the approaximation set, 
%where m is the number of objectives, n is the number of elements in the
%set. 
% 
%reference_point: The reference point of EHVI. A vector of length m. 
% 
% 
%
% Author: Guang Zhao <lemon15203@gmail.com>
% date: Oct 2017


% sort the front_array in one dimension
[~, I] = sort(front_array(:, 1));
front_array = front_array(I, :);

eh_integral = 0;
Phi = @(s) 1/2*(1+erf(s/sqrt(2)));
bigpsi = @(a, b, mu, sigma) sigma*normpdf((b-mu)/sigma)+(a-mu)*Phi((b-mu)/sigma);
[~, ncols]= size(front_array);

% odometer is used to exhaust every possible combination of [y1, y2, ...]
odometer_idx = ones(1, ncols);
odometer_bound = zeros(1, ncols);
odometer_value = cell(1, ncols);

% lower_corner upper_corner of each cell
lower_corner = zeros(1, ncols);
upper_corner = zeros(1, ncols);

for m =1:ncols
    odometer_value{m} = unique([front_array(:, m); -Inf; reference_point(m); Inf]);
    odometer_bound(m) = length(odometer_value{m})-1;
end

digit_idx = ncols;
while digit_idx
    for mm = 1:ncols
        lower_corner(mm) = odometer_value{mm}(odometer_idx(mm));
        upper_corner(mm) = odometer_value{mm}(odometer_idx(mm)+1);
    end
    
    %looking for active cell, whose lower bound dominates Yp, upper bound dominates r
    if point_dominat_set(front_array, lower_corner) &&...
            all(upper_corner<=reference_point) && any(upper_corner < reference_point)
            
            %compute contribution of each active cell
            delta_proder = zeros(1, ncols);
            fun = cell(1, ncols);
            for j = 1:ncols
                fun{j} = @(x) normcdf(x, mu(j), sigma(j));
            end
%             phi_proder = zeros(1, ncols);
%             v_corner = VCornerFind(ncols, odometer_idx, odometer_value, lower_corner, front_array, reference_point);
%             for j = 1:ncols
%                 delta_proder(j) = bigpsi(v_corner(j), upper_corner(j), mu(j), sigma(j))-...
%                     bigpsi(v_corner(j), lower_corner(j), mu(j), sigma(j));
%                 phi_proder(j) = Phi((upper_corner(j)-mu(j))/sigma(j))-Phi((lower_corner(j)-mu(j))/sigma(j));
%             end
%             S_minus = HyperSetCal(front_array, v_corner);
            
            for j = 1:ncols
                delta_proder(j) = integral(fun{j}, lower_corner(j), upper_corner(j));
            end
            
            delta = prod(delta_proder);
%             fun = @(x, y, z) normcdf(x, mu(1), sigma(1)).*normcdf(y, mu(2), sigma(2))...
%                 .*normcdf(z, mu(3), sigma(3));
%             xmin = lower_corner(1);
%             xmax = upper_corner(1);
%             ymin = lower_corner(2);
%             ymax = upper_corner(2);
%             zmin = lower_corner(3);
%             zmax = upper_corner(3);
%             
%             delta2 = integral3(fun,xmin,xmax,ymin,ymax,zmin,zmax);
%             
%             delta = prod(delta_proder)-S_minus*prod(phi_proder);
            if (delta<0)
                            
                hold on                
                scatter3(front_array(:,1), front_array(:,2), front_array(:,3), 'o')
                text(front_array(:,1), front_array(:,2), front_array(:,3), 'A')
                for k = 1:size(front_array, 2)
                    plot3([front_array(k, 1),reference_point(1)], [front_array(k, 2), front_array(k, 2)] , [front_array(k, 3), front_array(k, 3)])
                    plot3([front_array(k, 1),front_array(k, 1)], [front_array(k, 2), reference_point(2)] , [front_array(k, 3), front_array(k, 3)])
                    plot3([front_array(k, 1),front_array(k, 1)], [front_array(k, 2), front_array(k, 2)] , [front_array(k, 3), reference_point(3)])
                end
                scatter3(v_corner(:,1), v_corner(:,2), v_corner(:,3), '+')
                text(v_corner(:,1), v_corner(:,2), v_corner(:,3), 'v')
                scatter3(upper_corner(:,1), upper_corner(:,2), upper_corner(:,3),'*')
                text(upper_corner(:,1), upper_corner(:,2), upper_corner(:,3),'u')
                scatter3(0, 0, 0,'*')
                text(0, 0, 0,'l')
            


            end
            eh_integral = eh_integral+delta;%sum the contribution of all active cell
    end
    digit_idx = ncols;
    odometer_idx(digit_idx)  =  odometer_idx(digit_idx)+1;
    while odometer_idx(digit_idx) > odometer_bound(digit_idx)
        odometer_idx(digit_idx) = 1;
        digit_idx = digit_idx-1;
        if digit_idx == 0
            break
        end
        odometer_idx(digit_idx) = odometer_idx(digit_idx)+1;
    end
end

end

%%

function v_corner = VCornerFind(ncols, odometer_idx, odometer_value, lower_corner, front_array, reference_point)
    %Problem is how to define and find v, three dimension is quite
    %different with 2d.
    
    %C is a cell, delta is the contribution of C. S_minus is the region
    %dominated by every point in C
    v_corner = zeros(1, ncols);
    for m = 1:ncols
        temp_idx = odometer_idx(m);
        temp_point = lower_corner;

        while set_not_dominate_or_equal_to_point(front_array, reference_point, temp_point)
            temp_idx = temp_idx+1;
            temp_point(m) = odometer_value{m}(temp_idx);
        end
            
        
%         while point_dominate_set(front_array, temp_point)
%             temp_idx = temp_idx+1;
%             temp_point(m) = odometer_value{m}(temp_idx);
%         end
        v_corner(m) = odometer_value{m}(temp_idx);
    end
end



function set_volumn = HyperSetCal(S, reference_point)%return the hypervolumn of set S w.r.t. reference_point
    %based on algorithm in http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5766730
    %with objective of minimum 
    
    set_volumn = 0;
    
    for n = 1:size(S, 1)
        if all(S(n, :)< reference_point)
            set_volumn = set_volumn + ExclHvCal(S, n); 
        end
    end
    
    
    function exclhv = ExclHvCal(S, k)
        %S is a set, k is the idx of elements in S, the function returns
        %exclusive hypervolume of S(k) - S(k+1:end)
        limit_set = LimitSet(S, k);
        [limit_set_front, ~] = Find_pareto_front_multi2(limit_set);
        exclhv = InclHvCal(S(k, :)) - HyperSetCal(limit_set_front, reference_point);
    end
    
    function inclhv = InclHvCal(y)
        inclhv = prod(reference_point - y);
    end
    
    function limit_S = LimitSet(S, k)
%         S_minus = setdiff(S, y);
        limit_S = zeros(size(S, 1)-k, size(S, 2));
        for nn = k+1:size(S, 1)
            for mm = 1:size(S, 2)
                limit_S(nn-k, mm) = max(S(k, mm), S(nn, mm));
            end
        end
    end
end

function flag = set_not_dominate_or_equal_to_point(S, reference_point, y)
    flag = 1;
    
    for m = 1:length(y)
        if y(m) > reference_point(m)
            error('v doesnt dominate reference point')
        elseif y(m) == reference_point(m)
            flag = 0; % S and reference union already dominate y
            break;
        end
    end
    
    for n = 1:size(S, 1)
        if all(S(n, :) <= y)
            flag = 0;
            break;
        end
    end
end

% function flag = point_dominate_set(S, y)
% % if y dominates any elements in S, flag = 1
%     flag = 0;
%     for n = 1:size(S, 1)
%         if all(y <= S(n, :)) && any(y < S(n, :))
%             flag = 1;
%             break;
%         end
%     end
% end


function flag = point_dominat_set(S, y)
% if y is dominated by any elements in S, then flag = 0, otherwise flag = 1

flag = 1;

for n = 1:size(S, 1)
    if all(y >= S(n, :))% y is dominated or equal to m 
        flag = 0;
        break;
    end
end

end