clear
xi=[-7,-4,-3,0,3,4,7];
n=length(xi)
p0=ones(1,n)
pneg1= zeros(1,n)

t=1

w=ones(1,n);

for i=1:n
    for j=1:n
        if i==j
          w(i) = w(i)
        else
        w(i) = prod(w(i)/abs(xi(i)-xi(j)))
        end 
    end 
end

a=zeros(1,n);
b=zeros(1,n-1);
p=zeros(n-1,n);

%%%initial conditions 
a(1)= (xi.*p0.*w * p0')/(p0.*w * p0');
p(1,:)= (xi - a(1));
a(2)= ((xi.*p(1,:)) * p(1,:)'.*w)/(p(1,:)*p(1,:)'.*w);
b(1)= sqrt((p(1,:).*w * p(1,:)')/((p0.*w * p0')));
p(2,:)= ((xi - a(2)).*p(1,:)-b(1)^2*p0);

for i=3:n

    a(i) = (xi.*p(i-1,:).*w * p(i-1,:)') / (p(i-1,:).*w * p(i-1,:)');
    b(i-1) = sqrt(((p(i-1,:).*w*(p(i-1,:)')))/(p(i-2,:).*w*(p(i-2,:)')));
    p(i,:)=(xi-a(i)).*p(i-1,:)-b(i-1)^2*p(i-2,:);

end
%%constructing our jacobi matrix/hamiltonian

J= diag(a) + diag(b,-1) +diag(b,+1)
l=0.002
t=0:l:n;

for j=1:((n/l)+1)
    y(j)=U(J,t(j),n);
    x(j)=U2(J,t(j),n);
end

%%%ploting our graphs
hold on
plot(t, y,t,x)

annotation("line", [0.4789 0.4808], [0.9236 0.1111], "LineStyle", "--", "LineWidth", 2)
annotation("textbox", [0.4836 0.4307 0.1233 0.11], "String", "PST", "FontName", "Times New Roman", "FontSize", 36, "EdgeColor", "none")
 
grid on
legend(["|\langle e^{-iJt}e_0, e_N \rangle|", "|\langle e^{-iJt}e_0, e_0 \rangle|"], "Color", "none", "EdgeColor", "none", "FontSize", 24, "Position", [0.6769 0.2870 0.2519, 0.1703])
hAxes = findobj(gcf,"Type","axes")
hAxes.FontSize = 20
langleeiJte0e0rangle = findobj(gcf,"DisplayName","|\langle e^{-iJt}e_0, e_0 \rangle|")
langleeiJte0e0rangle.LineWidth = 2
langleeiJteNe0rangle = findobj(gcf,"DisplayName","|\langle e^{-iJt}e_0, e_N \rangle|")
langleeiJteNe0rangle.LineWidth = 2
hLegend = findobj(gcf,"Type","legend")
hLegend.EdgeColor = "none"
hLegend.Color = "none"
hLegend.FontSize = 24
hLineshape = findall(gcf,"Type","lineshape")
hLineshape.LineWidth = 2
hLineshape.LineStyle = "--"
hTextboxshape = findall(gcf,"Type","textboxshape")
hTextboxshape.FontName = "Times New Roman"
hTextboxshape.FontSize = 36
hTextboxshape.EdgeColor = "none"

%%%%%%Defining our functions
function fun= U(J1,S,n);

Ex= (expm(-i*J1*S));
fun= abs(Ex(1,n));

end

function fun2= U2(J2,S,n);


Ex2= (expm(-i*J2*S));
fun2=abs(Ex2(1,1));

end

