X=[ones(length(educ),1),educ,exper,smsa,black,south];
Y=lwage;
Beta=inv(X'*X)*X'*Y;
error=Y-X*Beta;
sigma=sum((error-mean(error)).^2)/(length(Y)-6);
covB=inv(X'*X)*sigma;
sigmaVar=2/(length(Y)-6)*sigma^2;
BVar=diag(covB);
VAR=[BVar;sigmaVar];
SigmaV=diag(VAR);
%theta=ones(7,1);
theta=[Beta;sigma];
A = normpdf(Y-X*theta(1:6),0,sqrt(theta(7)));
%L=prod(A);
L=sum(log(A));
i=1;
acc=0;
T=zeros(length(theta),10000);
r=zeros(length(theta),10000);
SigmaV_sc=SigmaV*0.1;
while i<10001;
    theta_new=theta+mvnrnd(zeros(7,1),SigmaV_sc)';
 while theta_new(7)<=0
      theta_new=theta+mvnrnd(zeros(7,1),SigmaV_sc)';
 end
A_new = normpdf(Y-X*theta_new(1:6),0,sqrt(theta_new(7)));
L_new=sum(log(A_new));
x=rand;
if exp(L_new-L)>=x
    acc=acc+1;
    L=L_new;
    theta=theta_new;
else
end
T(:,i)=theta;
i=i+1;
end
acc/10000
figure(1)
subplot(2,4,1)
hist(T(1,:),50)
title('beta 0')
subplot(2,4,2)
hist(T(2,:),50)
title('beta educ')
subplot(2,4,3)
hist(T(3,:),50)
title('beta exper')
subplot(2,4,4)
hist(T(4,:),50)
title('beta smsa')
subplot(2,4,5)
hist(T(5,:),50)
title('beta black')
subplot(2,4,6)
hist(T(6,:),50)
title('beta south')
subplot(2,4,7)
hist(T(7,:),50)
title('sigma')



mean(T,2);
X1 = norminv(0.975);
X2 = norminv(0.025);
CI=[0.035,0.085];
m=0.06;
S=(CI(2)-m)/X1;
theta=[Beta;sigma];
A = normpdf(Y-X*theta(1:6),0,sqrt(theta(7)));
L=sum(log(A))+log(normpdf(theta(2)-m,0,S));
i=1;
acc=0;
T=zeros(length(theta),10000);
SigmaV_sc=SigmaV*0.1;
while i<10001;
    theta_new=theta+mvnrnd(zeros(7,1),SigmaV_sc)';
 while theta_new(7)<=0
      theta_new=theta+mvnrnd(zeros(7,1),SigmaV_sc)';
 end
A_new = normpdf(Y-X*theta_new(1:6),0,sqrt(theta_new(7)));
L_new=sum(log(A_new))+log(normpdf(theta(2)-m,0,S));
x=rand;
if exp(L_new-L)>=x
    acc=acc+1;
    L=L_new;
    theta=theta_new;
else
end
T(:,i)=theta;
i=i+1;
end
acc/10000
figure(2)
subplot(2,4,1)
hist(T(1,:),50)
title('beta 0')
subplot(2,4,2)
hist(T(2,:),50)
title('beta educ')
subplot(2,4,3)
hist(T(3,:),50)
title('beta exper')
subplot(2,4,4)
hist(T(4,:),50)
title('beta smsa')
subplot(2,4,5)
hist(T(5,:),50)
title('beta black')
subplot(2,4,6)
hist(T(6,:),50)
title('beta south')
subplot(2,4,7)
hist(T(7,:),50)
title('sigma')





