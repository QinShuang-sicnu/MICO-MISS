%MIT License
%Copyright (c) [2021] [Shuang Qin]
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.
function x = MISS_SOCP(X,U,t,bias) %MISS_SOCP
% X is the sensor coordinates vector, e.g., [-20,-20;20,-20;20,20;-20,20]
% for 4 sensors.
% t is the measurement of each sensor
% bias is pre-set, here we use half of the value used in measurement
% generation
% U is Z or PZ in proposition 2
% The code numbers in the note are corresponding to the paper, u can check
% when needed
xr = [];
txr = [];
for i =1:size(U, 1)
    l = find(U(i,:) ==  1);
    m = find(U(i,:) == -1);
    if t(i) >= 0
        xr = [xr,[m;l]];
        txr = [txr,t(i)];
    else
        xr = [xr,[l;m]]; 
        txr = [txr,-t(i)];
    end
end
xc = unique(xr);
cvx_begin quiet
variable x(2) 
variable d(length(xc))
variable lam(length(xr))
minimize (norm(lam)+0.01*norm(d))
subject to
for m =1:length(xr)
    i = find(xc==xr(2,m));
    j = find(xc==xr(1,m));
    r = txr(m);
    si = X(xr(2,m),:)';
    sj = X(xr(1,m),:)';
    ca = r^2 + norm(sj)^2 - norm(si)^2  + 2*(si-sj)'*x;
    ca  + bias^2 + 2*bias*(r+d(j)) + 2*r*d(j) <= lam(m);%(25b)
    -ca - (bias^2 - 2*bias*(r+d(j)) + 2*r*d(j)) <= lam(m);%(25c)
%     -ca + r^2 + d(j)^2 <= lam(m);%(CASE 2)
    norm(x-sj) <= d(j);%(26)
end
cvx_end
end