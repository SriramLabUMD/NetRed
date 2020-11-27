function [ReducedModel] = NetRed(S,fluxes,mets,protectedMets,options)

% NetRed: An algorithm for metabolic network reduction.
%   This program is intended for use in reducing the results of Flux
%   Balance Analysis (FBA) simulations. For details of the algorithm,
%   please refer to the associated publication:
%   
%   Lugar, D., Mack, S., & Sriram, G. (2020). NetRed, an algorithm to reduce
%   genome-scale metabolic networks and facilitate the analysis of flux
%   predictions. Metabolic Engineering, doi: 10.1016/j.ymben.2020.11.003
%
%   If used, please cite the publication listed above.
%
% Inputs: NetRed(S,fluxes,mets,protectedMets,options)
%   S - The stoichiometric matrix
%   fluxes - (1) If analyzing a single flux solution: Enter a single flux
%                (column) vector solution corresponding to the given
%                stoichiometric matrix, with the fluxes in the same order
%                as their corresponding columns in S.
%            (2) If comparing multiple flux solutions: Enter a matrix whose
%                columns are composed of the flux solution vectors being
%                analyzed.
%   mets - A complete list of metabolite names in the same order as their
%          corresponding rows in S.
%   protectedMets - A list of metabolite names desired to remain in the
%          network by the user after network reduction.
%   options - User specified options may be included using the following
%          fields:
%          'filename' - Name of the file to which the results are written
%          'threshold' - Minimum value at which a flux is considered
%                   active (default 1E-8)
%          'lumpIdenticalRxns' - Flag (1 or 0) to combine parallel
%                   reactions into one after network reduction (default 0)
%          'verbose' - Flag (1 or 0) to print algorithmic details while
%                   reducing the network (default 1)
%          'BM' - Flag (1 or 0) to resolve degenerate biomass equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default variables
threshold = 10^-8;
lumpIdenticalRxns = 0;
filename = [];
verbose = 1;

% Handle options
if exist('options','var')
    if (isfield(options,'filename'))
        filename = options.filename;
    end
    if (isfield(options,'threshold'))
        threshold = options.threshold;
    end
    if (isfield(options,'lumpIdenticalRxns'))
        lumpIdenticalRxns = options.lumpIdenticalRxns;
    end
    if (isfield(options,'verbose'))
        verbose = options.verbose;
    end
    if (isfield(options,'BM'))
        BM = options.BM;
    end
end

% S is converted to a sparse matrix to decrease computation time.
S = sparse(S);
rIdx = (1:size(S,2))';

% Remove zero fluxes
zFlux = fluxes == 0; % zero fluxes
zFluxAll = all(zFlux,2); % fluxes that are zero in all flux
S(:,zFluxAll) = [];
fluxes(zFluxAll,:) = [];
rIdx(zFluxAll) = [];
zFlux(zFluxAll,:) = [];

% rows and columns are appended to the stoichiometric matrix corresponding
% to the dummy metabolites "MetIn" and "MetOut". All entry and exit fluxes
% must feed into and out of these dummy metabolites.
nMets = size(S,1);
mcol = zeros(nMets + 2,2);
mcol(nMets+1,1) = 1;
mcol(nMets+2,2) = -1;
mcol = sparse(mcol);

nRxns = size(S,2);
mrow = zeros(2,nRxns);
for j = 1:nRxns
    r = any(S(:,j) ~= 0,2);
    if sum(r) == 1
        if S(r,j) >= 0
            mrow(1,j) = -1*norm(S(r,j));
        end
        if S(r,j) <= 0
            mrow(2,j) = 1*norm(S(r,j));
        end
    end
end
S = [S;mrow];
S = [mcol,S];
fluxes(3:nRxns+2,:) = fluxes;
fluxes(1:2,:) = 0;
zFlux(3:nRxns+2,:) = zFlux;
zFlux(1:2,:) = 0;
rIdx = [0;0;rIdx];

% Reactions with negative flux values are negated such that their flux
% values are now positive.
nFluxes = size(fluxes,2);
Sset = cell(nFluxes,1);
[Sset{:}] = deal(S);
neg = fluxes < 0;
for i = 1:nFluxes
    Sset{i}(:,neg(:,i)) = -Sset{i}(:,neg(:,i));
    Sset{i}(:,zFlux(:,i)) = [];
end

protectedMets = unique(protectedMets);

% Finds the index of protected metabolites in the mets array.
pmets = find(ismember(mets,protectedMets)); %updated for speed -SM

% umets is the vector of indices of unprotected metabolites.
umets = 1:nMets;
umets(pmets) = [];
if size(umets,1) == 1
    umets = umets';
end

% Dummy metabolites MetIn and MetOut are appended to the metabolite list.
mets = [mets;'MetIn';'MetOut'];

if size(pmets,1) == 1
    pmets = pmets';
end
pmets = [pmets;[nMets+1;nMets+2]];

[TAll,umets,pmets,mustRemain] = findT(Sset,pmets,umets,mets,neg,zFlux,verbose);
Tfull = horzcat(TAll{:,1});
[~,iTa,iTc] = unique(Tfull','rows','stable');
T = Tfull(:,iTa);

bmKeepMets = [];
if exist('BM','var')
    BMnew = (rIdx == BM);
    for iter = 1:100
        lumpBM = find(T(BMnew,:));
        
        %Normalize lumped BM stoich so vBM_old = vBM_new
        T(:,lumpBM) = T(:,lumpBM)/diag(T(BMnew,lumpBM));
        
        nLumpBM = length(lumpBM);
        if nLumpBM > 1
            fprintf('Resolving degenerate biomass equations...\n')
            Tdiff = zeros(nRxns+2,1);
            for i = 2:nLumpBM
                Tdiff = Tdiff + T(:,lumpBM(1))-T(:,lumpBM(i));
            end
            Rdiff = find(abs(Tdiff)>1E-12);
            if isempty(Rdiff)
                T(:,lumpBM(2:nLumpBM)) = [];
                lumpBM = lumpBM(1);
                break
            end
            diffS = S(:,Rdiff);
            Tbm = T(Rdiff,lumpBM);
            mDiff = sum(diffS*Tbm,2);
            prodMet = find(mDiff>0);
            
            bmMet = setdiff(prodMet,pmets);
            if isempty(bmMet)
                reactMet = find(mDiff<0);
                bmMet = setdiff(reactMet,pmets);
            end
            bmKeepMets = [bmKeepMets; bmMet]; % set of unprotected products
            pmets = union(pmets,bmMet);
            umets = (1:nMets+2)';
            umets(pmets) = [];
            
            %Generate new T matrix
            [TAll,umets,pmets,tmp] = findT(Sset,pmets,umets,mets,neg,zFlux,verbose);
            mustRemain = [mustRemain;tmp];
            
            Tfull = horzcat(TAll{:,1});
            [~,iTa,iTc] = unique(Tfull','rows','stable');
            T = Tfull(:,iTa);
        elseif nLumpBM == 0
            keyboard;
        else
            break
        end
    end
end

Tct = zeros(size(T,2),1);
for i = 1:size(T,2)
    Tct(i) = sum(iTc == i);
end

% TO DO
%Scale T to ensure minimal stoichiometries
% for i = 1:size(T,2)
%     nz = T(:,i)~=0;
%     plus = all(abs(T(nz,i))>=1);
%     if plus
%         scale(i) = min(abs(T(nz,i)));
%     else
%         scale(i) = max(abs(T(nz,i)));
%     end
% end
% scale(lumpBM) = 1; %lumped BM eqn already scaled
% T = T./scale;

% Sred is the reduced stoichiometric matrix.
Sred = S*T;

if any(any(abs(Sred(umets,:))>1E-11))
    warning('Model reduction failed: Unprotected metabolites not properly removed')
end
Sred(umets,:) = 0;

% Removes infeasible loops from T matrix (empty columns in Sred)
loops = ~any(Sred,1);
if any(loops)
    T(:,loops) = [];
    Tct(loops) = [];
    Sred = S*T;
    Sred(umets,:) = 0;
    if exist('BM','var')
        lumpBM = find(T(BMnew,:));
    end
end

hit = false(size(T,2),nFluxes);
for i = 1:nFluxes
    hit(:,i) = ismember(T',TAll{i}','rows');
end
if exist('BM','var')
    hit(lumpBM,:) = 1;
end

% The following finds the new flux vector, newfluxes, for Sred. This is
% done by finding the minimum of the set of fluxes divided by their
% respective coefficient in the matrix T.

uT = 1./T;
uT(~isfinite(uT)) = 0;
newfluxes = cell(1,nFluxes);
for i = 1:nFluxes
    fmat = diag(fluxes(:,i))*uT;
    fmat(fmat == 0) = NaN;
    newfluxes{i} = min(fmat,[],1)';
    newfluxes{i}(isnan(newfluxes{i})) = 0;
end
newfluxes = horzcat(newfluxes{:});
newfluxes(~hit) = 0;

if lumpIdenticalRxns == 1
    tmp = (1:size(Sred,2))';
    [Sred,newfluxes,rMapRed,rMatRed] = combineparallel(Sred,newfluxes);
    if exist('BM','var')
        tmp = rMatRed*tmp;
        lumpBM = find(tmp==lumpBM);
    end
end

% Check to see if Sred*newfluxes equals the zero vector. If not, the
% command window will display a warning that the reduction has failed.
w = Sred*newfluxes;
w = w(1:nMets,:);
sumw = sum(abs(w));

stat = 0;
if any(sumw > threshold)
    warning('Network Reduction FAILED, newfluxes does not lie in the nullspace of Sred')
    stat = 1;
end

% ReducedModel is the structure containing the results of the network
% reduction.
ReducedModel.S = S; % Original stoichiometric matrix
ReducedModel.Sred = Sred; % Reduced stoichiometric matrix
ReducedModel.T = T; % Transition matrix
ReducedModel.Tct = Tct; % Transition matrix
ReducedModel.fluxes = fluxes; % Original set of fluxes (for S)
ReducedModel.newfluxes = newfluxes; % New set of fluxes (for Sred)
ReducedModel.pmets = sort(pmets); % List of protected metabolites
ReducedModel.umets = umets; % List of unprotected metabolites
ReducedModel.mets = mets; % List of all metabolites
ReducedModel.keepMets = mustRemain; % List of unremovable metabolites
ReducedModel.bmKeepMets = bmKeepMets; % List of metabolites kept to prevent BM degeneracy
ReducedModel.rxns = (1:size(Sred,2))'; % Number of each reaction
if exist('BM','var')
    ReducedModel.lumpBM = lumpBM; % Index of lumped biomass equation
end
ReducedModel.stat = stat;
ReducedModel.sumW = sumw;
if lumpIdenticalRxns == 1
    % TO DO
    %     ReducedModel.rMapBase = rMapBase; % Matrix containing sets of parallel/antiparallel in base model
    %     ReducedModel.rMatBase = rMatBase; % Matrix mapping representative reactions from base model
    ReducedModel.rMapRed = rMapRed; % Matrix containing sets of parallel/antiparallel in base model
    ReducedModel.rMatRed = rMatRed; % Matrix mapping representative reactions from base model
end

% Write the results to an Excel file.
if ~isempty(filename)
    s = {' '};
    arrow = '->';
    plus = '+';
    Sred = full(Sred);
    
    for j = 1:size(Sred,2)
        
        reactants = [];
        products = [];
        [u,~] = find(Sred(:,j) < 0);
        [a,~] = find(Sred(:,j) > 0);
        
        for i = 1:length(u)
            reactants = strcat(reactants,num2str(abs(Sred(u(i),j))),s,mets(u(i)),s);
            if i ~= length(u)
                reactants = strcat(reactants,s,plus,s);
            end
        end
        
        for k = 1:length(a)
            products = strcat(products,num2str(Sred(a(k),j)),s,mets(a(k)),s);
            if k ~= length(a)
                products = strcat(products,s,plus,s);
            end
        end
        
        reaction{j,1} = strcat(reactants,s,arrow,s,products);
    end
    
    for i = 1:length(reaction)
        reactionChar{i,1} = char(reaction{i});
    end
    
    warning('off','MATLAB:xlswrite:AddSheet')
    xlswrite(filename,ReducedModel.rxns,'Reactions','A1');
    xlswrite(filename,reactionChar,'Reactions','B1');
    xlswrite(filename,ReducedModel.newfluxes,'Reactions','C1');
    xlswrite(filename,ReducedModel.mets,'Metabolites','A1');
    if isempty(ReducedModel.umets) == 0
        xlswrite(filename,ReducedModel.umets,'RemovedMetabolites','A1');
    end
end
end

% BELOW ARE THE VARIOUS SUBROUTINES USED WITHIN THIS SCRIPT

function [r,c] = antiparallel(A)
% ANTIPARALLEL Determines the location of antiparallel vectors in the column space
% of a matrix
An = zeros(size(A));
for i = 1:size(A,2)
    if sum(A(:,i).^2) ~= 0
        An(:,i) = A(:,i)/norm(A(:,i));
    end
end
R = An'*An;
R = triu(R);
R = R - eye(size(R,1),size(R,2));
[r,c] = find(R <= -0.99999999);

end
function [r,c] = parallel(A)
% PARALLEL Determines the location of parallel vectors in the column space
% of a matrix
An = zeros(size(A));
for i = 1:size(A,2)
    if sum(A(:,i).^2) ~= 0
        An(:,i) = A(:,i)/norm(A(:,i));
    end
end
R = An'*An;
R = triu(R);
R = R - eye(size(R,1),size(R,2));
[r,c] = find(R >= 0.99999999);

end
function [Sp,vp,rMap,rMat] = combineparallel(S,v)
% Combines identical and reverse reactions in the stoichiometric matrix
% Inputs are S (the stoichiometric matrix) and v (the vector of fluxes)
% Outputs are Sc and vc, the stoichiometric matrix and flux vector with
% identical and reverse reactions combined.

vp = v;

[parPair(:,1),parPair(:,2)] = parallel(S);
[revPair(:,1),revPair(:,2)] = antiparallel(S);

rMap = speye(size(S,2));
for i = 1:size(parPair,1)
    p = parPair(i,:);
    if any(rMap(p(1),:))
        rMap(p(2),:) = rMap(p(2),:) + rMap(p(1),:);
        rMap(p(1),:) = 0;
        f = norm(S(:,p(2)))/norm(S(:,p(1)));
        for j = 1:size(vp,2)
            vp(p(2),j) = vp(p(2),j) + vp(p(1),j)/f;
        end
        vp(p(1),:) = 0;
    end
end
for i = 1:size(revPair,1)
    a = revPair(i,:);
    if all(any(rMap(a,:),2))
        rMap(a(2),:) = rMap(a(2),:) - rMap(a(1),:);
        rMap(a(1),:) = 0;
        f = norm(S(:,a(2)))/norm(S(:,a(1)));
        for j = 1:size(vp,2)
            vp(a(2),j) = vp(a(2),j) - vp(a(1),j)/f;
        end
        vp(a(1),:) = 0;
    end
end
rMat = speye(size(rMap));
tf = ~any(rMap,2);
rMat(tf,:) = [];
rMap(tf,:) = [];
Sp = S*rMat';
vp = rMat*vp;
end

function [targetRxns,unchangedRxns,StargetRxns,SunchangedRxns,StargetRxnsUmets,colMap] = submatricesS(S,umets)
% SUBMATRICESS Finds submatrices of S
% Inputs are S and umets

% targetRxns are the reactions involving unprotected metabolites.
targetRxns = any(S(umets,:) ~= 0); %streamlined generation and switched to logical indexing -SM

% unchangedRxns are those reactions which remain the same in the original
% and reduced S matrices. They do not involve unprotected metabolites.
unchangedRxns = 1:size(S,2);
unchangedRxns(targetRxns) = [];

% colMap maps the S matrix to one which separates unchanged rxns (left)
% from target rxns (right).
colMap = speye(size(S,2),size(S,2));
colMap = [colMap(:,unchangedRxns),colMap(:,targetRxns)];

% StargetRxns is the submatrix of S containing only columns which involve
% unprotected metabolites. SunchangedRxns is the submatrix of S containing
% only columns which do not involve unprotected metabolites.
StargetRxns = S(:,targetRxns);
SunchangedRxns = S(:,unchangedRxns);

% StargetRxnsUmets is the submatrix of StargetRxns containing only the rows
% of StargetRxns involving unprotected metabolites.
StargetRxnsUmets = S(umets,targetRxns);
end

function [PermMat,score] = findscore1num(S)
% FINDSCORE finds the "score" of each unprotected metabolite
%   A score greater than or equal to 0 is needed for a metabolite to be a
%   candidate for removal from the network.

S = full(S);
Mall = sum(S' < 0);
Nall = sum(S' > 0);
nmall = Nall'.*Mall';
score = Nall' + Mall' - nmall;

rowpermmat = eye(size(S,1));
rowpermmat = [score,nmall,rowpermmat];
rowpermmat = sortrows(rowpermmat);
rowpermmat = flipud(rowpermmat);

upperPermMat = rowpermmat(rowpermmat(:,1) >= 0,:);
lowerPermMat = rowpermmat(rowpermmat(:,1) < 0,:);
lowerPermMat(:,[1,2]) = [];

upperPermMat(:,1) = [];
upperPermMat = sortrows(upperPermMat);
upperPermMat(:,1) = [];

PermMat = sparse([upperPermMat;lowerPermMat]);

end

function [T2,mustRemain] = genvec5(S,umets,mets,verbose)
% Determines generating vector basis for the nullspace of a matrix (vectors of the pointed convex polyhedral cone)

zMets = ~any(S,2);
S = S(~zMets,:);
zMets = umets(zMets);
umets = setdiff(umets,zMets);

E = findscore1num(S);

S = E*S;

umets = E*umets;

St = S'; % St is the transpose of S

T = [St,eye(size(St,1))]; % This initiates matrix T for determining the generating vectors

mustRemain = [];
idx = 0:size(St,2)-1;
for iter = 1:100
    nMRi = length(mustRemain);
    m = 0;
    mustRemain = [];
    mr = [];
    for id = 1:length(idx)
        j = idx(id);
        
        A = T(:,(size(St,2)+1:size(T,2))); % This is the right hand side of the matrix T
        
        P = T(:,j+1)*T(:,j+1)'; %This outer product of T(j,:) says which entries in this vector have opposite signs from each other
        P = triu(P) - diag(diag(P));
        [row,col] = find(P < 0);
        
        set = [row,col]; % F1 are the pairs of elements of T(j,:) which have opposite signs
        set = sort(set,2);
        
        for k = 1:size(set,1)
            % Find intersection of the location of zeros in the first and second set rows
            I = all(A(set(k,:),:) == 0,1);
            
            %Make sure that this intersection is not a subset of the locations of zeros in the other rows
            Z = (A == 0);
            ii = zeros(1,size(A,2));
            ii(I) = 1;
            
            subset = ii*(Z');
            
            if sum(subset) ~= sum(ii)
                v = 1:size(A,1);
                v(set(k,:)) = [];
                
                if any(~any(A(v,I),2))
                    set(k,:) = [0 0]; % all set intersections which are subsets of other row zero locations are deleted
                end
            end
        end
        
        set = set(any(set,2),:); % removes rows of zeros in set
        
        set = sort(set,2);
        
        Tnew1 = T(abs(T(:,j+1)) < 1E-12,:); % Tnew1 are the rows of T that have T(:,j+1) near 0; switched to logical indexing for speed
        %         Tnew1 = T(T(:,j+1) == 0,:);
        
        Tnew2 = zeros(size(set,1),size(T,2));
        
        for i = 1:size(set,1)
            Tnew2(i,:) = abs(T(set(i,2),j+1))*T(set(i,1),:) + abs(T(set(i,1),j+1))*T(set(i,2),:); % Tnew2 are the new rows created by adding rows with entries of opposite signs
        end
        
        if rank(Tnew2) < size(Tnew2,1)
            m = m+1;
            mr(m,1) = j;
            mustRemain(m,1) = umets(j+1);
        else
            T = [Tnew1;Tnew2]; % This is a vertical concatenation, completing the development of vector T for this iteration
        end
    end
    idx = mr;
    nMRf = length(mustRemain);
    if nMRi == nMRf || nMRf == 0
        break
    end
end

T2 = T(:,(size(St,2)+1:size(T,2)));
T2 = T2'; % T2 are the generating vectors

if isempty(mustRemain) == 0 && verbose
    disp('The following metabolites must remain in the reduced network...')
    disp(mets(mustRemain))
end
end

function [TAll,umets,pmets,mustRemain] = findT(Sset,pmets,umets,mets,neg,zFlux,verbose)
% Finds all unremovable metabolites and then generates T matrices for each
% flux vector

fprintf('Finding all unremovable metabolites...\n')
mustRemain = [];
nFluxes = size(Sset,1);
[keep,Ttarget] = deal(cell(nFluxes,1));

% This loop lists any "unremovable" metabolites as protected, designated
% "mustRemain" on each iteration.
for j = 1:1000
    fprintf('Iteration %d \n',j)
    for i = 1:nFluxes
        StargetRxnsUmets = Sset{i}(umets,any(Sset{i}(umets,:) ~= 0));
        [Ttarget{i},keep{i}] = genvec5(StargetRxnsUmets,umets,mets,verbose);
        
        if ~isempty(keep{i})
            mustRemain = [mustRemain;keep{i}];
            pmets = [pmets;keep{i}];
            umets = (1:size(mets,1))';
            umets(pmets) = [];
        end
    end
    keepFull = vertcat(keep{:});
    
    if isempty(keepFull)
        break
    end
end

fprintf('Generating transition matrices for each flux...\n')
TAll = cell(nFluxes,1);
for i = 1:nFluxes
    
    % All important submatrices of S are found.
    [targetRxns,unchangedRxns,~,~,~,colMap] = submatricesS(Sset{i},umets);
    
    Ttarget{i} = [zeros(length(unchangedRxns),size(Ttarget{i},2));Ttarget{i}];
    Tunchanged = [eye(length(unchangedRxns));zeros(sum(targetRxns),length(unchangedRxns))];
    
    % This completes the generation of T (the transition matrix).
    Ttmp = [Tunchanged,Ttarget{i}];
    
    % This organizes the rows of T such that S*T = Sred.
    Ttmp = colMap*Ttmp;
    TAll{i} = sparse(size(zFlux,1),size(Ttmp,2));
    TAll{i}(~zFlux(:,i),:) = Ttmp;
    
    % Flip direction of reverse reactions
    TAll{i}(neg(:,i),:) = TAll{i}(neg(:,i),:)*-1;
    
end
end
