function collisions(data)
%calculates number of collisions per molecule for each dataset - uses 
%simplified output
output = [];
for(n=[1:size(data,2)])
    name=string(data(n).opts.name);
    collisions = sum(data(n).collisionsFullRaw) / sum(data(n).xEnd);
    output = [output;name,collisions];
end
output
end