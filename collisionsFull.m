function collisionsFull(Data)
%calculates number of collisions per molecule for each dataset - uses full
%output
output = [];
output2 = [];
for(n=[1:size(Data,2)])
    collisions = sum(sum(Data(n).collisionsFull > 0));
    collisions2 = sum(sum(Data(n).collisionsFull > 1));
    particles = sum(sum(~isnan(Data(n).xFull)));
    fraction = collisions / particles;
    output = [output;fraction];
    fraction2 = collisions2 / particles;
    output2 = [output2;fraction2];
end
output
output2
end