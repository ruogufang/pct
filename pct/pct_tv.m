function TVx = pct_tv(x,method)

if nargin<2
    method = 'iso';
end

TV_eps = 0;

% Total Variation norm of x, x is a n by n matrix
grad_x = [Grad1(x) Grad2(x) Grad3(x)];
if strcmp(method,'iso')
    pt_sqsum = sum(grad_x.*grad_x,2);
    if TV_eps == 0; TVx = sum(sqrt(pt_sqsum)); else TVx = sum(sqrt(pt_sqsum+TV_eps)); end
else
    TVx = norm(grad_x(:),1);
end

    function p=Grad1(u)
        % backward finite difference along dim 1
        p = cat(1,u(1,:,:),diff(u,1,1));
        p = p(:);
    end

    function q=Grad2(u)
        % backward finite difference along dim 2
        q = cat(2,u(:,1,:),diff(u,1,2));
        q = q(:);
    end

    function q=Grad3(u)
        % backward finite difference along dim 2
        q = cat(3,u(:,:,1),diff(u,1,3));
        q = q(:);
    end

end

