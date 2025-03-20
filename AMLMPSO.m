function [E_opt, A_opt, p_opt, S_opt, D_opt, X_reconstructed] = AMLMPSO(X, num_endmembers, alpha, beta, gamma)
    % Multi-Swarm PSO with Dimension-Division Strategy for AMLM Optimization
    % Swarm-wise optimization for each variable: E, A, p, S, D
    % Inputs:
    %   X - Hyperspectral image (M x N)
    %   num_endmembers - Number of endmembers (R)
    %   alpha, beta, gamma - Regularization parameters
    % Outputs:
    %   E_opt, A_opt, p_opt, S_opt, D_opt - Optimized variables
    %   X_reconstructed - Reconstructed hyperspectral image

    %% PSO Parameters
    numParticles = 20;  
    maxIter = 100;      
    w = 0.7;            
    c1 = 1.5;           
    c2 = 1.5;           

    %% Problem Dimensions
    [M, N] = size(X);
    R = num_endmembers;

    %% Initialize Global Best Solutions
    globalBest.E = VCA(X, R);
    globalBest.A = FCLS(X, globalBest.E);
    globalBest.p = rand(N, 1);      
    globalBest.S = rand(R, N);      
    globalBest.D = rand(M, R, N);   

    % Compute initial cost
    globalBest.cost = objectiveFunction(globalBest.E, globalBest.A, globalBest.p, globalBest.S, globalBest.D, X, alpha, beta, gamma);

    %% Initialize Five Separate Swarms
    swarms = struct();
    swarmVars = {'E', 'A', 'p', 'S', 'D'};
    
    for v = 1:length(swarmVars)
        varName = swarmVars{v};
        for i = 1:numParticles
            % Initialize variables with the global best solution
            swarm(i).(varName) = globalBest.(varName);
            
            % Initialize velocity
            swarm(i).V.(varName) = zeros(size(globalBest.(varName)));

            % Evaluate cost using fixed global best values for other variables
            swarm(i).cost = objectiveFunction(globalBest.E, globalBest.A, globalBest.p, globalBest.S, globalBest.D, X, alpha, beta, gamma);

            % Store personal best
            swarm(i).best = swarm(i);
        end
        swarms.(varName) = swarm;
    end

    %% Iterative Multi-Swarm PSO with Dimension-Division Strategy
    for iter = 1:maxIter
        for v = 1:length(swarmVars)
            varName = swarmVars{v};
            swarm = swarms.(varName);

            for i = 1:numParticles
                % Update particle velocity and position
                swarm(i) = updateParticle(swarm(i), globalBest, varName, w, c1, c2);

                % Evaluate new cost
                swarm(i).cost = objectiveFunction(globalBest.E, globalBest.A, globalBest.p, globalBest.S, globalBest.D, X, alpha, beta, gamma);

                % Update personal best
                if swarm(i).cost < swarm(i).best.cost
                    swarm(i).best = swarm(i);
                end
            end

            % **Dimension-Division Global Best Update**
            globalBest.(varName) = updateGlobalBest(globalBest.(varName), swarm, varName);

            % Store updated swarm
            swarms.(varName) = swarm;
        end
    end

    %% Extract Optimized Variables
    E_opt = globalBest.E;
    A_opt = globalBest.A;
    p_opt = globalBest.p;
    S_opt = globalBest.S;
    D_opt = globalBest.D;

    %% Reconstruct Hyperspectral Image
    P = ones(M, 1) * p_opt';
    Y = E_opt * (S_opt.*A_opt) + tensorProduct(D_opt, (S_opt.*A_opt));
    X_reconstructed = ((1 - P) .* Y) + (P .* Y) .* X;
    X_reconstructed = X_reconstructed./ sum(X_reconstructed,1);
end

%% **Dimension-Division Global Best Update Function**
function newGlobalBest = updateGlobalBest(globalBest, swarm, varName)
    % Initialize the new global best with the current values
    newGlobalBest = globalBest;
    
    % Number of particles
    numParticles = length(swarm);
    
    % Get the size of the variable
    dimSize = size(globalBest);
    
    % Preallocate storage based on variable dimensions
    if isvector(globalBest)  % For vector variables (e.g., p)
        allValues = zeros(length(globalBest), numParticles);
    elseif ndims(globalBest) == 2  % For matrix variables (E, A, S)
        allValues = zeros([dimSize, numParticles]);
    elseif ndims(globalBest) == 3  % For tensor D
        allValues = zeros([dimSize, numParticles]);
    else
        error('Unexpected variable dimensions.');
    end

    % Store costs for all particles
    allCosts = zeros(1, numParticles);

    % Fill in the values for each particle
    for i = 1:numParticles
        allCosts(i) = swarm(i).best.cost;  % Store the cost function value
        if isvector(globalBest)
            allValues(:,i) = swarm(i).best.(varName);
        elseif ndims(globalBest) == 2
            allValues(:,:,i) = swarm(i).best.(varName);
        elseif ndims(globalBest) == 3
            allValues(:,:,:,i) = swarm(i).best.(varName);
        end
    end
    
    % Find the index of the best particle
    [~, bestIdx] = min(allCosts);
    
    % Assign the best values per dimension
    if isvector(globalBest)
        newGlobalBest = allValues(:,bestIdx);
    elseif ndims(globalBest) == 2
        newGlobalBest = allValues(:,:,bestIdx);
    elseif ndims(globalBest) == 3
        newGlobalBest = allValues(:,:,:,bestIdx);
    end
end

%% Updated Cost Function Including γΩ(D)
function cost = objectiveFunction(E, A, p, S, D, X, alpha, beta, gamma)
    % Compute Derived Variables
    M = size(E,1);
    N = size(A,2);
    
    P = ones(M, 1) * p';
    Y = (E * (S .* A)) + tensorProduct(D, S .* A);  
    
    % Compute Reconstruction Error
    X_hat = ((1 - P) .* Y) ./ (1 - P .* Y);
    reconstructionError = norm(X - X_hat, 'fro')^2;  
    
    % Compute Total Variation Regularization
    TV_A = totalVariation(A);
    TV_S = totalVariation(S);

    % Compute Energy Constraint for D
    energy_D = sum(arrayfun(@(i) norm(D(:,:,i), 'fro')^2, 1:N));

    % Compute Final Cost Function
    cost = reconstructionError + alpha * TV_A + beta * TV_S + gamma * energy_D;
end

%% Particle Update Function
function particle = updateParticle(particle, globalBest, varName, w, c1, c2)
    % Update velocity for only the selected variable
    particle.V.(varName) = w * particle.V.(varName) + ...
                            c1 * rand() .* (particle.best.(varName) - particle.(varName)) + ...
                            c2 * rand() .* (globalBest.(varName) - particle.(varName));

    % Update position for only the selected variable
    particle.(varName) = particle.(varName) + particle.V.(varName);
end


%% Función para calcular la variación total (TV)
function tv = totalVariation(matrix)
    tv = sum(abs(diff(matrix, 1, 1)), 'all') + sum(abs(diff(matrix, 1, 2)), 'all');
end

%% Producto tensorial de D con S .* A
function result = tensorProduct(D, SA)
    [M, R, N] = size(D);
    result = zeros(M, N);
    for j = 1:N
        result(:, j) = D(:, :, j) * SA(:, j);
    end
end

%%
function Po = VCA(Y,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P,indices,SNRe]=VCA(Y,N)
%
% Vertex Component Analysis algorithm for endmembers estimation in multi/hyperspectral dataset
%  
%
% Inputs
%   Y --> Multi/hyperspectral dataset as 2D matrix (L x K).
%   N --> Number of endmembers to find.
%
% Outputs
%   P --> Matrix of endmembers (L x N).
%
% References
%   J. M. P. Nascimento and J. M. B. Dias, ?Vertex component analysis: A 
% fast algorithm to unmix hyperspectral data,? IEEE Transactions on 
% Geoscience and Remote Sensing, vol. 43, no. 4, apr 2005.
%
% DUCD February/2021
% IICO-FC-UASLP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization.
K = size(Y, 2);
L = size(Y, 1);

yMean = mean(Y, 2);
RZeroMean = Y - repmat(yMean, 1, K);
[Ud, ~, ~] = svds(RZeroMean*RZeroMean.'/K, N);
Rd = Ud.'*(RZeroMean);
P_R = sum(Y(:).^2)/K;
P_Rp = sum(Rd(:).^2)/K + yMean.'*yMean;
SNR = abs(10*log10( (P_Rp - (N/L)*P_R) / (P_R - P_Rp) ));

SNRth = 15 + 10*log(N) + 8;
if (SNR > SNRth) 
    d = N;
    [Ud, ~, ~] = svds((Y*Y.')/K, d);
    Yd = Ud.'*Y;
    u = mean(Yd, 2);
    M =  Yd ./ repmat( sum( Yd .* repmat(u,[1 K]) ) ,[d 1]);
else
    d = N-1;
    r_bar = mean(Y.').';
    Ud = pca(Y, d);
    %Ud = Ud(:, 1:d);
    R_zeroMean = Y - repmat(r_bar, 1, K);
    Yd = Ud.' * R_zeroMean;
     c = zeros(N, 1);
    for j=1:K
        c(j) = norm(Yd(:,j));
    end
    c = repmat(max(c), 1, K);
    M = [Yd; c];
end
e_u = zeros(N, 1);
e_u(N) = 1;
A = zeros(N, N);
% idg - Doesnt match.
A(:, 1) = e_u;
I = eye(N);
%k = zeros(K, 1);
for i=1:N
    w = rand(N, 1);
    % idg - Oppurtunity for speed up here.
    tmpNumerator =  (I-A*pinv(A))*w;
    %f = ((I - A*pinv(A))*w) /(norm( tmpNumerator ));
    f = tmpNumerator / norm(tmpNumerator);

    v = f.'*M;
    k = abs(v);
    [~, k] = max(k);
    A(:,i) = M(:,k);
    indices(i) = k;
end
if (SNR > SNRth)
    Po = Ud*Yd(:,indices);
else
    Po = Ud*Yd(:,indices) + repmat(r_bar, 1, N);
end
return;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function [U] = pca(X, d)
    N = size(X, 2);
    xMean = mean(X, 2);
    XZeroMean = X - repmat(xMean, 1, N);     
    [U,~,~] = svds((XZeroMean*XZeroMean.')/N, d);
return;
end

%%
function A = FCLS(X, E_init)
    % Fully Constrained Least Squares (FCLS) Unmixing
    % Solves for abundance matrix A given endmembers E_init and hyperspectral image X
    % Inputs:
    %   - E_init: Matrix of endmembers (M x R), where M is the number of spectral bands
    %   - X: Matrix of hyperspectral pixels (M x N), where N is the number of pixels
    % Output:
    %   - A: Abundance matrix (R x N), where R is the number of endmembers
    
    R = size(E_init,2); % Number of spectral bands and endmembers
    N = size(X,2);      % Number of pixels
    
    % Preallocate abundance matrix
    A = zeros(R, N);
    
    % Constraints: Non-negative and sum-to-one
    lb = zeros(R, 1); % Lower bound (non-negative abundances)
    Aeq = ones(1, R); % Equality constraint (sum-to-one)
    beq = 1;
    options = optimoptions('lsqlin', 'Algorithm', 'interior-point', 'Display', 'off');
    
    % Solve for each pixel
    for i = 1:N
        A(:, i) = lsqlin(E_init, X(:, i), [], [], Aeq, beq, lb, [], [], options);
    end
end

