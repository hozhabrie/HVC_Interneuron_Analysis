function array = treeClusts(outperm,bounds)
    for i = 1:length(bounds)
        [~,loc1] = find(outperm == bounds(i,1));
        [~,loc2] = find(outperm == bounds(i,2));
        array{i} = sort(outperm(loc1:loc2));
    end
end