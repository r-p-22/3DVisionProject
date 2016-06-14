filenameVectorsBeforeBA = 'width_vectors_before_BA.txt';
filenameVectorsAfterTrivialBA = 'width_vectors_0_0.txt';
filenameVectorsAfterAugmentedBA = 'width_vectors_10_100.txt';

[vectorsBeforeBA, ~] = importdata(filenameVectorsBeforeBA)
[vectorsAfterTrivialBA, ~] = importdata(filenameVectorsAfterTrivialBA)
[vectorsAfterAugmentedBA, ~] = importdata(filenameVectorsAfterAugmentedBA)

figure

scatter3(vectorsBeforeBA(:,1),vectorsBeforeBA(:,2),vectorsBeforeBA(:,3),'MarkerFaceColor','green', 'MarkerEdgeColor','black');
 
hold on;

scatter3(vectorsAfterTrivialBA(:,1), vectorsAfterTrivialBA(:,2), vectorsAfterTrivialBA(:,3),'MarkerFaceColor','blue', 'MarkerEdgeColor','black');

hold on;

scatter3(vectorsAfterAugmentedBA(:,1), vectorsAfterAugmentedBA(:,2), vectorsAfterAugmentedBA(:,3),'MarkerFaceColor','red', 'MarkerEdgeColor','black');

view(80,30);
axis equal;
