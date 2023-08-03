function J = imbin(I, bin, dim)
    if ~isscalar(dim) %% multidimension bin
        assert(numel(bin)==numel(dim))
        for n = 1:numel(dim)
            I = imbin(I, bin(n), dim(n));
        end
        J = I;
    else %% single dimension bin
        D = size(I);
        L = D(dim);
        L = L - mod(L, bin);
        %% cut part of image that cannot be binned
        S.type = '()';
        S.subs = {};
        for nD = 1:numel(D)
            if nD~=dim
                sliceidc = 1:D(nD);
            else
                sliceidc = 1:L;
            end
            S.subs = [S.subs, sliceidc];
        end
        I = subsref(I, S);
        %% bin by mean
        D2 = [D(1:dim-1), bin, L/bin, D(dim+1:end)];
        D3 = [D(1:dim-1), L/bin, D(dim+1:end)];
        I = reshape(I, D2);
        I = mean(I, dim);
        I = reshape(I, D3);
        J = I;
    end

