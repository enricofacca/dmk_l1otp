% Routine for substituting a point to the closest one in the Euclidean
% distance sense
function p=substitute_point_in_mesh(p,fixed_point);
% From old coordinates matrix p_old onw gets the new coordinates matrix
% p_new by substituting fixed_point to its closest point in p_old
% fixed_points must be given in horizontal (N\times 3 matrix)
distances = zeros(size(p,1),size(fixed_point,1));

for j=1:size(fixed_point,1)
    for i=1:size(p,1)
        distances(i,j) = norm(p(i,:)-fixed_point(j,:));
    end
    if sum(distances(:,j)==min(distances(:,j))) == 1
        p(distances(:,j)==min(distances(:,j)),:)
        p(distances(:,j)==min(distances(:,j)),:) = fixed_point(j,:);
    else
        disp('Error on the distances')
        break;
    end
end
