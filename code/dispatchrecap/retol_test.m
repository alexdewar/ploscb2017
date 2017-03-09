function retol_test

close all

[A,B,C,D,E,F,G,H] = deal(true(3));
disp('all')
A(1,:) = false;
B(1,:) = false;
retol(A,B);
retol(B,A);

disp('none')
C(1,:) = false;
D(2,:) = false;
retol(C,D);
retol(D,C);

disp('third')
E(1,:) = false;
F(:,1) = false;
retol(E,F);
retol(F,E);

disp('another')
G(1,:) = false;
H(1,1:2) = false;
retol(G,H);
retol(H,G);

dump2base(true)

function rol=retol(A,B)
if ~nargin
    A = true(2);
    A(1,:) = false;
    B = true(2);
    B(1,:) = false;
end

A = ~A;
B = ~B;

Av = A(:);
Bv = B(:);
overlap = Av & Bv;
Q = sum(overlap);
rol = 0.5 * (Q/sum(Av) + Q/sum(Bv))

figure
subplot(1,2,1)
imshow(~A)
subplot(1,2,2)
imshow(~B)

function rol=retol2(A,B)
if ~nargin
    A = true(2);
    A(1,:) = false;
    B = true(2);
    B(1,:) = false;
end

A = ~A;
B = ~B;

Av = A(:);
Bv = B(:);
overlap = Av & Bv;
Q = sum(overlap);
R = sum(~Av & Bv);

rol = Q/(Q+R)

figure
subplot(1,2,1)
imshow(~A)
subplot(1,2,2)
imshow(~B)
