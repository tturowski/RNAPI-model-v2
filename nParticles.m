function nParticles(data)
%calculates number of collisions per molecule for each dataset - uses 
%simplified output
output = [];
for(n=[1:size(data,2)])
    name=string(data(n).opts.name);
    nParticles = sum(data(n).nParticles);
    output = [output;name,nParticles];
end
output
end