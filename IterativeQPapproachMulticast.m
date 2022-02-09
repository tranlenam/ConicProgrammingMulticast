clear variables
clc
%********************************************************************
% problem specs
%********************************************************************nTx=8;
% Number of transmit antennas
nTx = 1000; 
% Number of users
nUsers = 100; 
% Maximum number of iterations
nIterations = 100;  
% Number of samples used in the warm-up rpcoder
nWarmups = 1000;
% Select the options for YALMIP
ops = sdpsettings('solver','mosek','verbose',0,'debug',0); % set the interal solver to be Mosek
w = sdpvar(nTx,1,'full','complex'); % complex beamformer to be computed
obj= w'*w; %the objective to be minimized, i.e., (6a)

% Generate channel realizations.
channel = sqrt(1/2)*(randn( nUsers, nTx)+1i*randn(nUsers, nTx));
% Warm-up procedure. We generate nWarmups random matrices and pick up  the
% one with smallest norm
wnextwarmup = randn(nTx,nWarmups) + 1i*randn(nTx,nWarmups);
wnextwarmup = wnextwarmup./repmat(min(abs(channel*wnextwarmup),[],1),nTx,1);
[~,idx] = min(diag(wnextwarmup'*wnextwarmup));
wnext = wnextwarmup(:,idx);
x = real(channel*w); % see the definition of (6c) in the paper
y = imag(channel*w); % see the definition of (6c) in the paper
% Set the intial point
x0 = real(channel*wnext);
y0 = imag(channel*wnext);
seqobj = zeros(nIterations,1); % the sequence of objectives
errtol = 1e-2;
%% Iterative procedure now starts 
for iIteration = 1:nIterations
    F = [(x0.^2 + 2*x0.*(x-x0) + y0.^2 + 2*y0.*(y-y0)) >= ones(nUsers,1)];
    diagnostics = solvesdp(F,obj,ops); % solve the problem
    if(diagnostics.problem == 0)
        x0 = double(real(channel*w)); 
        y0 = double(imag(channel*w));
        seqobj(iIteration) = (double(w'*w));
        % Break if decrease over 10 consecutive iterations is less than errtol 
        if((iIteration>10)&&(abs(seqobj(iIteration) - seqobj(iIteration-10)) < errtol))
            seqobj(min(iIteration+1,nIterations):end) = [];
            break;
        end
    else
        % print error if solver can't solve the approximate problem
        disp(yalmiperror(diagnostics.problem))
        seqobj = NaN;
        break;
    end
end
% Compute SNR of each user, should be larger 1 for all users
AchievedSNR = abs(channel*double(w)).^2;
disp(['The required power is ',num2str(10*log10(seqobj(end))),' dBw'])
%Plot the convergence of the objective
plot(1:length(seqobj),10*log10(seqobj))
xlabel('Iteration count, $n$','Fontsize',12,'Interpreter','latex') % add the label to the x-axis
ylabel('Required Power $20\log_{10}||\mathbf{w}||_{2}$','Fontsize',12,'Interpreter','latex') % add the label to the y-axis
%saveas(gcf, '../../results/ConvergencePlot.png')
