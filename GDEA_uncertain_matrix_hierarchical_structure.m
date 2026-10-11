% This code puts together the hierarchical matrices for an AHP problem, 
% where pairwise evaluations are made using intervals rather than crisp numbers, 
% and calculates the corresponding constraints to be added to DEA-AHP model
% After running this code and evaluating the alternatives on the bottom lever criteria
% you can run dea_ahp code to find the efficiency scores for all the alternatives
% For more details see the header of the dea_ahp.m file

% The first part of the code is asking you for the structure of the hierarchical composition of your problem.
% The code assumes that your hierarchical structure has at least 3 levels.
% The first hierarchical level is considered to have 1 element by default.
L       = input('\nHow many hierarchical levels in this problem? >> ');
C       = eye(L,1);
C(2)    = input('\nHow many criteria in level 2? >> ');
% Construct the hierarchical structure for the 3-level GDEA model.
    I = [0; 0];
    for i=2:L-1
        m   = 0;
        k   = sum(C(i,:));
        s   = sprintf('\n======================= Level %d ========================',i);
        disp(s);
        for j=1:k
            if ~isempty(find(I(i,:)==j))
                C(i+1,j)        = 1;
                I(i+1,j)        = sum(C(i+1,:));
                m = m+1;
            else
                str             = ['How many sub-criteria for the criterion ',num2str(j-m),' in level ',num2str(i),'? >> '];
                e               = input(str);
                if e == 0
                    C(i+1,j)    = 1;
                    I(i+1,j)    = sum(C(i+1,:));
                else
                    C(i+1,j)    = e;
                    I(i+1,j)    = 0;
                end
            end
        end
    end

    % Calculating matrix P that will help position all the submatrices within a big constraint matrix
    P           = C(2:L,:);
    for i=L:-1:3
        Ind         = [0 C(i-1, find(C(i-1,:)))];
        for j=1:size(Ind,2)-1
            P(i-2,j)    = sum(P(i-1,sum(Ind(1:j))+1:sum(Ind(1:j+1))));
        end
    end
    
    % The process below is designed to let user enter only the neccessary elements for each pariwise
    % matrix and for each hierarchical level
    C0              = C;
    C0(find(C0==1)) = 0;
    [i,j,v]         = find(C0);
    q               = [i j v];
    pm              = sortrows(q,[1 2]);
    nL              = zeros(L,1);
    for i=2:L
        nL(i)           = size(find(pm(:,1)==i),1);
    end
    
    disp(' ');
    disp('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp('  %%%  When entering pairwise matrix elements,  %%%');
    disp('  %%%    enter only the neccessary elements     %%%');
    disp('  %%%      from the upper triangular part,      %%%');
    disp('  %%%     enclose the elements in brackets,     %%%');
    disp('  %%%   and separate them by spaces or commas.  %%%');
    disp('  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    
    % initilazing the resulting AHP constraint matrix h:
    h               = [];
    
    % within the following for loop a section of the matrix h is calculated for each hierarchical level
    for i=2:L-1
        width           = sum(C(i,:));
        s               = sprintf('\n================= Hierarchical Level %d =================',i);
        disp(s);
        
        s = sprintf('Enter the pairwise matrix for Level %d:\n', i);
        disp(s);
        mi              = find(pm(:,1) == i);
        mim             = pm(mi,3);
    
        for j = 1:nL(i)
            s               = sprintf('Level %d, Pairwise Matrix %d (dimensions are %d x %d):',i,j,mim(j),mim(j));
            disp(s);
            a               = ones(mim(j),mim(j));
            w               = ones(mim(j),mim(j));
    
            % Prompting you to enter the elements for the current matrix at the current hierarchical level
            % Entering is done by rows; one row of the pairwise matrix at the time
            disp('First enter all the lower values:');
            for k=1:mim(j)-1
                str             = ['Row ',num2str(k),' : '];
                a(k,k+1:mim(j)) = input(str);
            end 
            disp('Now enter all the upper values for the same matrix:');
            for k=1:mim(j)-1
                str             = ['Row ',num2str(k),' : '];
                w(k,k+1:mim(j)) = input(str);
            end 
    
            % Given a pairwise upper triangular matrix a with lower triangular part being all 1
            % this part of the program calculates the corresponding constraint matrix based on all the
            % existing relationships between the two criteria that can be established using
            % the maximum number of two steps 
            n                       = size(a,1);
            bl                      = triu(a)+tril(1./a')-eye(size(a));
            bu                      = triu(w)+tril(1./w')-eye(size(w));
            Pairwise_Matrix_Lower   = bl
            Pairwise_Matrix_Upper   = bu
            k                       = 0;
            aq                      = [];
            for ii=1:n-1
                for jj=ii+1:n
                    k                   = k+2;
                    apmin               = [(bl(ii,:).*bu(:,jj)') (bl(ii,:).*bl(:,jj)')];
                    apmax               = [(bu(ii,:).*bl(:,jj)') (bu(ii,:).*bu(:,jj)')];
                    a1                  = min(apmin);
                    a2                  = -max(apmax);
                    aq(k-1:k,[ii jj])   = [-1 a1; 1 a2];
                end
            end
            Lower_and_Upper_Limits  = aq
            qq{j}           = aq;
            sizeaq(j,:)     =size(aq);
        end
    
        % Initialising the constraint matrix Q{i} for level i and populating it by the small
        % constraint matricies calculated above using each pairwise matrix from the same hierarchical
        % level in turn
        iC0             = find(C0(i,:));
        Q{i}            = zeros(sum(sizeaq(:,1)),width);
        iq              = 1;
        iaq             = 1;
        jq              = 0;
        jaq             = 0;
        z               = 0;
        for j=1:width
            if ~isempty(find(iC0==j))
                z                                                                   = z+1;
                Q{i}(iq:sum(sizeaq(1:iaq,1)),j+jq-jaq:j+jq-jaq+sizeaq(jaq+1,2)-1)   = qq{z};
                iq                                                                  = iq + sizeaq(iaq,1);
                iaq                                                                 = iaq+1;
                jq                                                                  = sum(sizeaq(1:jaq+1,2));
                jaq                                                                 = jaq+1;
            end
        end
        % Q{i} constraint matrix corresponding to the i-th hierarchical level is now calculated
    
        % Need to appropriately widen the matrix Q{i} to fit the width of the widest matrix Q{i}
        if i~=L
            y           = [];
            for j=1:size(Q{i},2)
                y           = [y nncopy(Q{i}(:,j),1,P(i,j))];
            end
        else
            y           = Q{i};
        end
        % The widened matrix Q{i} is now stored in matrix y and that matrix in the next sted i
        % appended to the previous version of matrix h:
        h           = [h; y];
    end

    % Finally adding one more column to the matrix h (to account for the efficiency variable in dea_g
    h           = [h zeros(size(h,1),1)];

% And saving the matrix h into the file called ahph.
save ahph h;
disp('#################################################################');
disp('## Constraints for DEA_AHP are now calculated. If you have     ##');
disp('## also evaluated all alternatives with respect to all the     ##');
disp('## bottom level criteria you can now run dea_ahp from Matlab   ##');
disp('#################################################################');
disp(sprintf(' '));
