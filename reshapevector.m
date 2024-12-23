function [v] = reshapevector(v)

global problemID nballs N rankY

if ( problemID == 2 )

    xa = [];
    xb = [];
    xa(:,1) = v(1:2*nballs);
    xa = reshape(xa, nballs, [])';
    xb.uv = mat2cell(xa, 2, ones(1, nballs))';
    xb.s = v(2*nballs+1:3*nballs);
    xb.r = v(end);
    v = [];
    v = xb; 

elseif ( problemID == 3 )

    A = reshape(v(1:4), [2, 2]);
    b = reshape(v(5:6), [2, 1]);

    v = [];
    v.A = A;
    v.b = b;

end