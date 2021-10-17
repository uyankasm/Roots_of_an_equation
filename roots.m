clear;
clear all;

e0=((1/(36*pi))*10^-9);...Define epsilon
iterasyon=1;
lower=-3;
upper=10;
...define equation
syms 'x'
denklem(x)=(1/(4*pi*e0))*((13*(x+7)/(abs(x+7))^3)+(9*(x+4)/(abs(x+4))^3)+(6*(x-11)/(abs(x-11))^3)+(3*(x-14)/(abs(x-14))^3));
...define derivative
dfx(x)=81000000000/abs(x + 4)^3 + 117000000000/abs(x + 7)^3 + 54000000000/abs(x - 11)^3 + 27000000000/abs(x - 14)^3 - (27000000000*sign(x + 4)*(9*x + 36))/abs(x + 4)^4 - (27000000000*sign(x - 14)*(3*x - 42))/abs(x - 14)^4 - (27000000000*sign(x - 11)*(6*x - 66))/abs(x - 11)^4 - (27000000000*sign(x + 7)*(13*x + 91))/abs(x + 7)^4;
   
......................................................................

...BİSECTİON METHOD...
    
mid=(lower+upper)/2;
error=abs((upper-lower)/(2^iterasyon));
   fprintf("\n------------------THE BISECTION METHOD------------------------\n");
   fprintf("I\t\t\tERROR\t\t\t\t\t\tX\t\t\t\t\tF(X)\n");
   
while error> 10^(-10)
    iterasyon=iterasyon+1;
    if denklem(mid)*denklem(upper) <0
        lower=mid;
    else
        upper=mid;
    end
    mid=(lower+upper)/2;
    error=abs((upper-lower)/(2^iterasyon));
    grafikb(iterasyon)=error;
 
    fprintf("%d\t\t%e\t\t %.20f \t %e\n",iterasyon,error,mid,denklem(mid));
end
   

...................................................................... 
    
...NEWTON METHOD...
iterasyon1=0;
lower=-3;
upper=10;
x1=(lower+upper)/2;
x0=x1+1;

 fprintf("\n--------------------THE NEWTON METHOD---------------------\n");
 fprintf("I\t\t\tERROR\t\t\t\t\t\tX\t\t\t\t\tF(X)\n");
 while abs(x1-x0)>10^-10
    iterasyon1=iterasyon1+1;
     x0=x1;
     x1=x0-(denklem(x0)/dfx(x0));
     grafikn(iterasyon1)=abs(x1-x0);
    fprintf("%d\t\t%e\t\t %.20f \t %e\n",iterasyon1,abs(x1-x0),x1,denklem(x1));
 end
 
.....................................................................
    
% %  ...SECANT METHOD
%       
% lower=-3;
% upper=10;
% p0=(upper+lower)/2;
% p1=p0+1;
% p=0;
% while abs(p-p1)< 10^-10
%   f0=denklem(p0);
%   f1=denklem(p1);
%     p=p1-((f1*(p1-p0))/(f1-f0));
%     p0=p1;
%     p1=p;
% end
%     fprintf("The root is with The secant methot");
%     format long
%     double(p)
%     double(denklem(p))

......................................................................
    

...Secant method v1.2 ...  
p0=(upper+lower)/2;
p1=p0+1;

fprintf("\n--------------------THE SECANT METHOD---------------------\n");
 fprintf("I\t\t\tERROR\t\t\t\t\t\tX\t\t\t\t\tF(X)\n");
 for i=1:1000
 f0=denklem(p0);
 f1=denklem(p1);
y=p1-(f1*(p1-p0)/(f1-f0));
err=abs(y-p1);
grafiks(i)=err;
fprintf("%d\t\t%e\t\t %.20f \t %e\n",i,err,y,denklem(y));
if err<10^-10
break
end
p0=p1;
p1=y;
 end
 

    xb=1:iterasyon;
    xn=1:iterasyon1;
    xs=1:i;
    semilogy(xb,grafikb,xn,grafikn,xs,grafiks);
    xlabel('iteration')
    ylabel('error')
    
    legend('Bisection','Newton','Secant')



