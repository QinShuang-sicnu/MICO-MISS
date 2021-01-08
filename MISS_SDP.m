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
function x = MISS_SDP(X,U,t,bias) %MISS_SDP
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
cvx_begin sdp quiet
variable y(2+length(xc))
variable Y(2+length(xc),2+length(xc))  hermitian 
variable lam(length(xr))
minimize (sum(lam))
for m =1:length(xr)
    i = find(xc==xr(2,m));
    j = find(xc==xr(1,m));
    r = txr(m);
    si = X(xr(2,m),:)';
    sj = X(xr(1,m),:)';
    c1 = r^2+norm(sj)^2-norm(si)^2+bias^2+2*bias*r;
    c2 = r^2+norm(sj)^2-norm(si)^2+bias^2-2*bias*r; 
    u1 = [2*(si-sj);zeros(length(xc),1)];
    u1(2+j) = 2*(bias+r);
    u2 = [2*(si-sj);zeros(length(xc),1)];
    u2(2+j) = 2*(-bias+r);
    Q = [eye(2),zeros(2,length(xc));zeros(length(xc),2+length(xc))];
    Q(2+j,2+j) = -1;
    subject to
    c1^2 + trace(u1*u1'*Y) + 2*c1*u1'*y <= lam(m)%(29a)
    c2^2 + trace(u2*u2'*Y) + 2*c2*u2'*y <= lam(m)%(29b)
    trace(Q*Y) - 2*sj'*y(1:2) + norm(sj)^2 == 0 %(30)
    norm(y(1:2)-sj) <= y(2+j);%(26)
end
[Y,y;y',1] >= 0;%(31)
cvx_end
x = y(1:2);
end