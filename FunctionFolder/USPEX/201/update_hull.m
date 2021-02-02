function hull =update_hull( hull,E,X,mu_min, mu_max, IND)
X10 = min(hull(:,1));
X20 = max(hull(:,1));
imin = 0;
imax = 0;
np = size(hull,1);
good = 0;
%update the already known hull point, from left to right
for i = 1:np
    if abs(X-hull(i,1))<0.001  %stable point
       if E<hull(i,2)
          good = 1;
       else
          return
       end
    elseif X<hull(i,1)
        if hull(i,1) < X20 + 0.001
           X2 = i;
           X20 = hull(i,1);
        end
        imax = imax + 1;
    elseif X>hull(i,1)
        if hull(i,1) > X10 - 0.001
           X1 = i;
           X10 = hull(i,1);
        end
        imin = imin + 1;
    end
end

if ~good
    E0= hull(X1,2)+  (X-X10)*(hull(X1,2)-hull(X2,2))/(X10-X20);
    if E>E0 %
       return
    else
       good = 1;
    end
end
todo = zeros(1,np);

%from right side
for i = imin:-1:1
  mu = (E-hull(i,2)) / (X-hull(i,1));
  if i==1
     if mu<mu_min
        hull(i,2) = E-mu_min*(X-hull(i,1));
     end
  else
     mu1 = (hull(i,2)-hull(i-1,2))/(hull(i,1)-hull(i-1,1));
     if mu<mu1
        todo(i)=1;  %to update        
     else
        break;
     end
  end
end
%from left side
if imin+imax < np
    imin = imin+1;
    todo(imin)=1;
end

for i =(imin+1):1:(imin+imax)
  mu = (hull(i,2)-E) / (hull(i,1)-X);
  if i==np
     if mu>mu_max
        hull(i,2) = E-mu_max*(X-hull(i,1));
     end
  else
     mu1 = (hull(i,2)-hull(i+1,2))/(hull(i,1)-hull(i+1,1));
     if mu>mu1
        todo(i)=1;  %to update
     else
        break;
     end
  end
end

%to update
count = 1;
for i=1:imin
    if todo(i)<1
       hull1(count,:)=hull(i,:);
       count = count + 1;
    end
end

hull1(count,:)=[X,E,IND];
for i=(imin+1):(imin+imax)
    if todo(i)<1
       count = count +1;
       hull1(count,:)=hull(i,:);       
    end
end
hull = hull1;
