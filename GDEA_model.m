%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Before running this code all the alternatives have to be evaluated on %
%  bottom level criteria and the scores have to be stored in one matrix  % 
%  called BOLECREVA. Each row of BOLECREVA corresponds to an alternative,%
%  and each column to one bottom level criterion. Store the matrix into  %
%  the file using the command:                                           %
%                 save filename BOLECREVA                                %
%  where "filename" is any name you choose. Alternatively, you can run:  %
%  writeinmatrix code (type writeinmatrix in the MATLAB command prompt), %
%  which will allow you to enter the matrix values and it will           %
%  automatically be saved in the file named 'blmatrix'.                  %
%  Also, before running this code, you need to calculate the AHP         %
%  constraint matrix using ahp_pairwise command or ahp_uncertain. Use    %
%  ahp_uncertain in case the entries in pairwise matrices are intervals. % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% asking for the filename where BOLECREVA matrix is stored:
BLC = input('\nWhat is the file name where BOLECREVA is stored?\n>> ','s');
% loading that file and reading the matrix BOLECREVA into memory:
load(BLC);
% loading ahph file where the AHP constraint matrix was stored after 
% executing ahp_pairwise or ahp_uncertain command. 
% ahph contains the constraint matrix whose name is h.
load ahph;
hr                  = size(h,1);
% finding the size of the BOLECREVA matrix; n will be the number of rows, 
% and io is the number of columns:
[n,io]              = size(BOLECREVA);
% taking logarithm of all the elements of the BOLECREVA matrix 
% to bring it into the form required by DEA-G:
M                   = log(BOLECREVA);
% setting the objective function coefficients
f                   = zeros(io+1,1);
f(io+1,1)           = -1;
% setting the LP constraint matrices: A and Aeq
A                   = zeros(n+hr,io+1);
A(1:n,1:io)         = M;
A(1:n,io+1)         = ones(n,1);
A(n+1:end,:)        = h;
Aeq                 = ones(1,io+1);
Aeq(1,io+1)         = 0;
% setting the RHS vectors: b and beq
b                   = zeros(n+hr,1);
beq                 = 1;
% setting the lower bounds on the decision variables:
vlb                 = zeros(io+1,1);
vlb(io+1,1)         = -Inf;
% setting the upper bounds on the decision variables:
vub                 = Inf*ones(io+1,1);
% saving anotgher copy of the matrix A into B:
B                   = A;
% initiating the vector lp_geom where the results will be stored
lp_geom             = zeros(n,1);
% for each alternative making neccessary changes to the constraint 
% matrix A and finding the optimal solution
for j=1:n,
    A                   = B;
    Ak                  = ones(n,1)*A(j,1:io);
    A(1:n,1:io)         = A(1:n,1:io)-Ak;
    x                   = linprog(f,A,b,Aeq,beq,vlb,vub);
    lp_geom(j,1)        = f'*x;
    X(:,j)              = x;
end
% the weights each alternative has chosen are stored in matrix X, 
% which is being saved in file weights:
save weights X;
% displaying the weights calculated (one column for each alternative):
weights_for_each_alternative    = X(1:end-1,:)
% and the best scores obtained for each alternative; 
% 1 is the highest score; 0 is the lowest:
dea_ahp_results                 = exp(-lp_geom)'