function backtracking(Data)
%calculates fraction of bactracked RNAP - uses full
%output
output = [];
for(n=[1:size(Data,2)])
    backtracked = sum(sum(Data(n).backtrackingFull <= -1));
    particles = sum(sum(~isnan(Data(n).xFull)));
    fraction = backtracked / particles;
    output = [output;fraction];
end
output
end