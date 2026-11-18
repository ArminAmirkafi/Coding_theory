function msg_cap = CA_SCL(rq, n, nL, F, info_pos, K, A, crcg)

    N = 2^n;
    f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b)); % minsum
    g = @(a,b,c) b + (1-2*c).*a;

    LLR = zeros(nL, n+1, N); 
    u_cap = zeros(nL, n+1, N); 
    
    PML = inf * ones(nL,1);
    PML(1) = 0;
    node = 0;
    depth = 0; 
    done = 0;    
    ns = zeros(1, 2*N-1); 
    u_cap(:, n+1, node+1) = 0;
    
    LLR(:, 1, :) = repmat(rq, nL, 1, 1);
    
    while done == 0
        if depth == n
            DM = LLR(:, n+1, node+1);
            DM = reshape(DM, nL, 1); 
            
            if any(F == node+1)
                u_cap(:, n+1, node+1) = 0; 
                PML = PML + abs(DM).*(DM < 0);
            else
                dec = DM < 0; 
                PM2 = [PML; PML + abs(DM)];
                [PML, pos] = mink(PM2, nL); 
                
                pos1 = pos > nL; 
                pos(pos1) = pos(pos1) - nL;
                
                dec = dec(pos); 
                dec(pos1) = 1 - dec(pos1); 
                
                LLR = LLR(pos, :, :);
                u_cap = u_cap(pos, :, :); 
                u_cap(:, n+1, node+1) = dec;
            end
            
            if node == N-1
                done = 1;
            else
                node = floor(node/2);
                depth = depth - 1;
            end
            
        else 
            npos = (2^depth-1) + node + 1;
            
            if ns(npos) == 0 % left child
                temp = 2^(n-depth);
                Ln = LLR(:, depth+1, temp*node+1:temp*node+temp);
                Ln = reshape(Ln, nL, temp);
                
                a = Ln(:, 1:temp/2);
                b = Ln(:, temp/2+1:end);
                temp = temp/2; 
                node = node*2;
                depth = depth+1;
                
                LLR(:, depth+1, temp*node+1:temp*node+temp) = f(a,b);
                ns(npos) = 1;
            
            elseif ns(npos) == 1 % right child
                temp = 2^(n-depth);
                Ln = LLR(:, depth+1, temp*node+1:temp*node+temp);
                Ln = reshape(Ln, nL, temp);
                
                a = Ln(:, 1:temp/2);
                b = Ln(:, temp/2+1:end);
                 
                left_node = 2*node;
                left_depth = depth+1;
                left_temp = temp/2;
                
                u_cappp = u_cap(:, left_depth+1, left_temp*left_node+1:left_temp*(left_node+1));
                u_cappp = reshape(u_cappp, nL, left_temp);
                
                node = 2*node + 1;
                depth = depth + 1; 
                temp = temp/2;
                
                LLR(:, depth+1, temp*node+1:temp*(node+1)) = g(a,b,u_cappp);
                ns(npos) = 2;
                
            elseif ns(npos) == 2 % return to parent
                temp = 2^(n-depth);
                lnode = 2*node;
                rnode = 2*node + 1;
                child_depth = depth + 1; 
                child_temp = temp/2;
                
                u_capl = u_cap(:, child_depth+1, child_temp*lnode+1:child_temp*(lnode+1));
                u_capl = reshape(u_capl, nL, child_temp);
                
                u_capr = u_cap(:, child_depth+1, child_temp*rnode+1:child_temp*(rnode+1));
                u_capr = reshape(u_capr, nL, child_temp);
                
                u_cap(:, depth+1, temp*node+1:temp*(node+1)) = [mod(u_capl+u_capr,2), u_capr];
                
                node = floor(node/2);
                depth = depth - 1;
            end
        end
    end
    
    msg_cap1 = u_cap(:, n+1, info_pos);
    msg_cap1 = reshape(msg_cap1, nL, K);
    
    %path selection
    cout = 1; 
    for cl = 1:nL
        [~, r1] = gfdeconv(fliplr(msg_cap1(cl,:)), crcg);
        if isequal(r1, 0)
            cout = cl;
            break;
        end
    end
    msg_cap = msg_cap1(cout, 1:A);
end