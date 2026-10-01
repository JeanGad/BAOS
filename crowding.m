function cd = crowding(F)
% CROWDING  NSGA-II crowding distance of one front.
[n, M] = size(F); cd = zeros(n,1); if n <= 2, cd(:) = inf; return; end
for m = 1:M
  [fs, o] = sort(F(:,m)); rg = fs(end) - fs(1); if rg == 0, rg = 1; end
  cd(o(1)) = inf; cd(o(end)) = inf;
  cd(o(2:end-1)) = cd(o(2:end-1)) + (fs(3:end) - fs(1:end-2))/rg;
end
end
