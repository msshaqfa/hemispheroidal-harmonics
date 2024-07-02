% Bijection between the unit disk and the Northern hemispheroid with radii a,a,c.
function P = spheroidal_projection(p,a,c)
    if size(p, 2) == 1
      p = [real(p), imag(p)];
    end
    u = p(:,1);
    v = p(:,2);
    if size(p,2) < 3
      z = 1 + u.^2 + v.^2;
      P = [2*a*u./z, 2*a*v./z, -c*(-1+u.^2+v.^2)./z];
      
      P(isnan(z)|(~isfinite(z)),1) = 0;
      P(isnan(z)|(~isfinite(z)),2) = 0;
      P(isnan(z)|(~isfinite(z)),3) = -c;
        
    else
      z = p(:,3);
      P = [(u/a)./(1+z/c), (v/a)./(1+z/c)];
      P(isnan(P)) = Inf;
    end
end
