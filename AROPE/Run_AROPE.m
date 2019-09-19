% Refer to AROPE.m for details
A=load('datasets/drugsim1network.txt'); % change the path to the network you want to embed
order = [1,2,3,-1];
weights = cell(4,1);
weights{1} = 1;
weights{2} = [1,0.1];
weights{3} = [1,0.1,0.01];
weights{4} = 0.001;
[psim4U_cell,V_cell] = AROPE(A,50,order,weights); % you can change the embedding dimension here
save('dsim1U_cell.mat','dsim1U_cell');