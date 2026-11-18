function u_hat = SC_decoder(l,K,Q)

N = length(l);
n = log2(N);

frozen = Q(1:end-K);
is_frozen = false(1, N);
is_frozen(frozen) = true;

L = zeros(n+1, N);
L(1,:) = l; 
U_hat = zeros(n+1, N);
u_hat=U_hat(n+1,:);
 

f = @(a,b) (1-2*(a<0)).*(1-2*(b<0)).*min(abs(a),abs(b)); 
g = @(a,b,c_val) b+(1-2*c_val).*a; 

node = 0;
depth = 0;
done = 0;
ns = zeros(1, 2*N-1);
while done==0
    temp=N/2^depth; 
    npos=2^depth+node;
    startIdx=temp*node+1;
    endIdx=temp*node+temp;
    
    if depth ==n
        idx=node+1;
        if is_frozen(idx)
            u_hat(idx)=0;
        else
            u_hat(idx) = double(L(depth+1,idx)<0); 
            
            if L(depth+1,idx) < 0
            u_hat(idx) = 1;
            else
            u_hat(idx) = 0;
            
            end
        end
        U_hat(depth+1,idx)=u_hat(idx);
        
        if depth==0
            done=1;

        else
            node=floor(node/2);
            depth=depth-1;
        end
    else
        L0=L(depth+1,startIdx:endIdx);
        a=L0(1:temp/2);
        b=L0(temp/2+1:end);
        
        if ns(npos) == 0
            LLR_l=f(a,b);
            sub_start=startIdx;
            sub_end=startIdx+temp/2-1;
            L(depth+2,sub_start:sub_end)=LLR_l;
            node = 2*node;
            depth = depth+1;
            ns(npos) = 1;
  
        elseif ns(npos) == 1
            sub_start_left = startIdx;
            sub_end_left = startIdx+temp/2-1;
            c_val = U_hat(depth+2, sub_start_left:sub_end_left);
            
            LLR_r = g(a,b,c_val);
            sub_start_right = startIdx+temp/2;
            sub_end_right = endIdx;
            L(depth+2, sub_start_right:sub_end_right) = LLR_r;
            ns(npos)=2;
            node=2*node+1;
            depth=depth+1;
            
        elseif ns(npos) == 2
            sub_start_left = startIdx;
            sub_end_left = startIdx+temp/2-1;
            sub_start_right = startIdx+temp/2;
            sub_end_right = endIdx;
            
            left_part=U_hat(depth+2, sub_start_left:sub_end_left);
            right_part=U_hat(depth+2, sub_start_right:sub_end_right);
            U_hat(depth+1,startIdx:endIdx)=[mod(left_part+right_part,2), right_part];
            
            if depth==0
                done=1;
            else
                node=floor(node/2);
                depth=depth-1;
            end
        end
    end
end
u_hat=U_hat(n+1,:);

end
