function [A,b,c,lbounds,ubounds,FEASIBLE] = preprocess...
         (A,b,c,lbounds,ubounds,BIG)
% PREPROCESS  - Preprocessing input data.
% Usage:  [A,b,c,lbounds,ubounds] = preprocess(A,b,c,lbounds,ubounds,BIG)

% Yin Zhang, last updated April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global NNZA
global b_orig c_orig
global data_changed
global Lbounds_non0
global Ubounds_exist nub
global Fixed_exist ifix infx xfix
global Zrcols_exist izrcol inzcol xzrcol
global Sgtons_exist isolved insolved xsolved

fprintf('Preprocessing ...\n');
[m, n] = size(A); FEASIBLE = 1;

if any(lbounds > ubounds)
   fprintf('\nPreprocessor: Lower bound exceeds upper bound\n');
   fprintf(1,'%c',7);  % ring a bell
   FEASIBLE = 0; return;
end;

if ~issparse(A) A = sparse(A); end;
b = sparse(b); c = sparse(c);
lbounds = sparse(lbounds);
b_orig = b; c_orig = c;
data_changed = 0;



%----- find upper bounds -----
nub = 0; iubounds = ubounds < BIG - 1;
Ubounds_exist = full(any(iubounds));
if (Ubounds_exist)
   ubounds = sparse(iubounds.*(ubounds-lbounds)); 
   nub = nnz(ubounds);
end

[m, n] = size(A); NNZA = nnz(A);
fprintf(' (m=%i, n=%i)\n',m,n);

% Foram eliminadas todas as operacoes do preprocessamento
% Ficaram sparce, limites superiores
