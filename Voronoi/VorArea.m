
function A = VorArea(V,C)

A = zeros(length(C),1);

for i = 1:length(C)
    A(i,1) = polyarea(V(C{i},1),V(C{i},2));
end

% A = A';

end