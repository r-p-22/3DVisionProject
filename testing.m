points = csvread('group18good.csv');
cands = csvread('candidates18good.csv');



s=[];

for i=1:size(cands,1)
    
    D = cands';
    D(:,i) =[];
    D = [D -cands(i,:)'];
    [U,S,V] = svd(D);

    sol = V(:,end)/V(end,end);

    s = [s sol];
end

smod = mod(s,1);
isnotint = smod > 0.05 & smod < 0.95;
find(sum(isnotint)==0)

bvecs = csvread('finalbasis18good.csv');
bounds = csvread('boundaries18.csv');


figure,
scatter3(points(:,1),points(:,2),points(:,3)),hold on;
for i=1:2
    plot3(points(1,1)+[0;bvecs(i,1)],...
        points(1,2)+[0;bvecs(i,2)],points(1,3)+[0;bvecs(i,3)]);hold on;
end
scatter3(bounds(:,1),bounds(:,2),bounds(:,3),'rs'),hold on;


figure,
for i=1:size(cands,1)
    plot3([0;cands(i,1)],[0;cands(i,2)],[0;cands(i,3)]);hold on;
end


