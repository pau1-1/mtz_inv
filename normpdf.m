function pdf = normpdf (x, m, s)

if (nargin ~= 1 && nargin ~= 3)
    error('normpdf: you must give one or three arguments');
end

if (nargin == 1)
    m = 0;
    s = 1;
end

if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
        error ('normpdf: x, m and s must be of common size or scalars');
    end
end

sz = size (x);
pdf = zeros (sz);

if (isscalar (m) && isscalar (s))
    if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
        pdf = NaN * ones (sz);
    else
        pdf = stdnormal_pdf ((x - m) ./ s) ./ s;
    end
else
    k = find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf));
    if (any (k))
        pdf(k) = NaN;
    end

    k = find (~isinf (m) & ~isnan (m) & (s >= 0) & (s < Inf));
    if (any (k))
        pdf(k) = stdnormal_pdf ((x(k) - m(k)) ./ s(k)) ./ s(k);
    end
end

pdf((s == 0) & (x == m)) = Inf;
pdf((s == 0) & ((x < m) | (x > m))) = 0;

end
