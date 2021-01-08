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
function y = MICO(X,U,t,bias) %MTCO
% X is the sensor coordinates vector, e.g., [-20,-20;20,-20;20,20;-20,20]
% for 4 sensors.
% t is the measurement of each sensor
% bias is pre-set, here we use half of the value used in measurement
% generation
% The code numbers in the note are corresponding to the paper, u can check
% when needed
n = length(X);%number of sensors
U=[U,U]; % U corresponding to equation (4) in the paper
cvx_begin sdp quiet
variable H(2*n,2*n)  hermitian
variable Y(2,2)  hermitian 
variable h(2*n)
variable y(2)
minimize ( real ( trace(t*t' - 2*U*h*t' + U*H*U') ...
    + 10^(-4)*sum(sum(H)) ) )
subject to 
for i =1:n
    h(i) >= 0;%(10a) 
    h(i) <= bias;%(10b)
    real(H(i,i)) >= h(i)^2;%(13)
    real(H(i,i)) <= bias^2;%(11a)
    for j = 1:n
        if(i==j)
            continue
        else
            real(H(i,j)) >= 0; 
            real(H(i,j)) <= bias^2;%(11b)
        end
    end
end
for i =n+1:2*n
    h(i) >= norm(X(i-n,:)'-y);%(12)
    real(H(i,i)) == trace(Y) - 2*X(i-n,:)*y + X(i-n,:)*X(i-n,:)';%(9b)
    for j = n+1:2*n
        if(i==j)
            continue
        else
        real(H(i,j)) >= abs(trace([eye(2,2),-X(j-n,:)';...
            -X(i-n,:),X(i-n,:)*X(j-n,:)']...
            *[Y,y;y',1]));%(16)
        end
    end
end
[Y,y;y',1] >= 0;%(9c)
[H,h;h',1] >= 0;%(9c)
cvx_end
end