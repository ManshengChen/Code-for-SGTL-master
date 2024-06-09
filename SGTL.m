function [S0,obj,resofeachiter]=SGTL(H,c,y,alpha,beta,omega)
V = size(H,3);
N = size(H,2); %sample number

for v=1:V
    S{v} = eye(N);
    Z{v} = eye(N);
    Y1{v} = zeros(N,N);
    Y2{v} = zeros(N,N);
    A{v} = eye(N,N);
end

Isconverg = 0;epson = 1e-4;
iter = 0;
mu = 10e-5; max_mu = 10e10; pho_mu = 2;
rho = 0.0001; max_rho = 10e12; pho_rho = 2;
tic;

foreachiter = 1;
resofeachiter=[];
sX = [N, N, V];

while(Isconverg == 0)
    %fprintf('----processing iter %d--------\n', iter+1);
    for v=1:V
        % update Z^k and S^k
        B1 = Z{v} - Y1{v}/mu;
        B2 = A{v} - Y2{v}/rho;
        B12 = (-H(:,:,v) + mu * B1 + rho * B2)/(rho + mu);
        B = (B12 + B12')/2;
        S{v} = project_fantope(B,c);
        Z{v} = prox_l1((S{v}+Y1{v}/mu),alpha/mu);
        Y1{v} =  Y1{v} + mu*(S{v} - Z{v});
    end
    
    
    % update G
    S_tensor = cat(3, S{:,:});
    Y2_tensor = cat(3, Y2{:,:});
    s = S_tensor(:);
    y2 = Y2_tensor(:);
    [a, ~] = wshrinkObj(s+1/rho*y2,beta/rho,sX,0,3,omega);
    A_tensor = reshape(a, sX); 
    
    
    % update W
    Y2_tensor = Y2_tensor + rho*(S_tensor - A_tensor);
    
    
    %% coverge condition
    Isconverg = 1;
    maxvalue = 0;
    for v=1:V
        A{v} = A_tensor(:,:,v);
        Y2{v} = Y2_tensor(:,:,v);
        if (norm(S{v}-A{v},inf)>epson)
            history.norm_S_A = norm(S{v}-A{v},inf);
            maxvalue = max(maxvalue,history.norm_S_A);
            Isconverg = 0;
        end
        if (norm(S{v}-Z{v},inf)>epson)
            history.norm_S_Z = norm(S{v}-Z{v},inf);
            maxvalue = max(maxvalue,history.norm_S_Z);
            Isconverg = 0;
        end
    end
    history.objval(iter+1)=maxvalue;
    if (iter>30)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
    rho = min(rho*pho_rho, max_rho);
    
    if foreachiter
        S0 = zeros(N);
        for v=1:V
            S0 = S0 + abs(Z{v})+abs(Z{v}');
        end
        results=clustering8(abs(S0)+abs(S0'),c, y);
        resofeachiter = [resofeachiter;results];
    else
        resofeachiter = [];
    end
end
S0 = zeros(N);
for v=1:V
    S0 = S0 + abs(Z{v})+abs(Z{v}');
end
S0 = S0 -diag(diag(S0));
obj =history.objval;
