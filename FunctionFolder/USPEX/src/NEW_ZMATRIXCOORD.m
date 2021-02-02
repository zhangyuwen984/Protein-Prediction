function Coords= NEW_ZMATRIXCOORD(ZMATRIX,format)

if(size(ZMATRIX,1)==1)
Coords=ZMATRIX;
else
TEEMPI = ZMATRIX(1,:);
ZMATRIX(1,:)=zeros(1,3);
Coords = zeros(0,3);
Coords(1,:)=[0,0,0];
Coords(2,:)=[ZMATRIX(2,1),0,0];
if size(ZMATRIX,1)>2.5
    for ind = 3:size(ZMATRIX,1)
        if ind==3
            RelatedC=[ZMATRIX(format(ind,1),1),0,0;ZMATRIX(format(ind,2),1),0,0;0,0,1];
        else

            form(3) = format(ind,3);
            form(2) = format(ind,2);
            form(1) = format(ind,1);
           
            RelatedC= Coords(form,:);
        end

        r=ZMATRIX(ind,1);
        theta=ZMATRIX(ind,2);
        phi=ZMATRIX(ind,3);
        [Coords(ind,1),Coords(ind,2),Coords(ind,3)]= conStructCOO(RelatedC,r,theta,phi);

        if ind==3
            %two rotations
            for outer=1:2
                if outer==1
                    
                    rotVec = [0,1,0];
                    torsang=ZMATRIX(2,2)-pi/2;
                    
                else
                    rotVec=[0,0,1];
                    torsang=ZMATRIX(2,3);

                end

                for inder=1:3

                    endVector = [cos(torsang)+rotVec(1)^2*(1-cos(torsang))  , rotVec(1)*rotVec(2)*(1-cos(torsang)) - rotVec(3)*sin(torsang), rotVec(1)*rotVec(3)*(1-cos(torsang)) + rotVec(2)*sin(torsang);...
                        rotVec(2)*rotVec(1)*(1-cos(torsang)) + rotVec(3)*sin(torsang), cos(torsang)+rotVec(2)^2*(1-cos(torsang)) , rotVec(2)*rotVec(3)*(1-cos(torsang)) - rotVec(1)*sin(torsang);...
                        rotVec(3)*rotVec(1)*(1-cos(torsang)) - rotVec(2)*sin(torsang), rotVec(3)*rotVec(2)*(1-cos(torsang)) + rotVec(1)*sin(torsang), cos(torsang)+rotVec(3)^2*(1-cos(torsang)) ]*Coords(inder,:)';
                    Coords(inder,:)=endVector';
                end

            end
        end
    end
    
    
elseif size(ZMATRIX,1)==2
      for outer=1:2
                if outer==1
                    rotVec = [0,1,0];
                    torsang=ZMATRIX(2,2)-pi/2;

                else

                    rotVec=[0,0,1];
                    torsang=ZMATRIX(2,3);
                end


                for inder=1:2


                    endVector = [cos(torsang)+rotVec(1)^2*(1-cos(torsang))  , rotVec(1)*rotVec(2)*(1-cos(torsang)) - rotVec(3)*sin(torsang), rotVec(1)*rotVec(3)*(1-cos(torsang)) + rotVec(2)*sin(torsang);...
                        rotVec(2)*rotVec(1)*(1-cos(torsang)) + rotVec(3)*sin(torsang), cos(torsang)+rotVec(2)^2*(1-cos(torsang)) , rotVec(2)*rotVec(3)*(1-cos(torsang)) - rotVec(1)*sin(torsang);...
                        rotVec(3)*rotVec(1)*(1-cos(torsang)) - rotVec(2)*sin(torsang), rotVec(3)*rotVec(2)*(1-cos(torsang)) - rotVec(1)*sin(torsang), cos(torsang)+rotVec(3)^2*(1-cos(torsang)) ]*Coords(inder,:)';
                    Coords(inder,:)=endVector';
                end

       end
    
end

Coords(:,1)=Coords(:,1)+TEEMPI(1,1);
Coords(:,2)=Coords(:,2)+TEEMPI(1,2);
Coords(:,3)=Coords(:,3)+TEEMPI(1,3);
end

Coords=real(Coords);

function [x,y,z]= conStructCOO(RelatedC,r,theta,phi)

RelatedC = RelatedC';

xi = RelatedC(1,1);
yi = RelatedC(2,1);
zi = RelatedC(3,1);
xj = RelatedC(1,2);
yj = RelatedC(2,2);
zj = RelatedC(3,2);
xk = RelatedC(1,3);
yk = RelatedC(2,3);
zk = RelatedC(3,3);
%
%  Find unit vector along j->i vector
%
xji = xi - xj;
yji = yi - yj;
zji = zi - zj;
rji = xji*xji + yji*yji + zji*zji;
rji = sqrt(rji);
xji = xji/rji;
yji = yji/rji;
zji = zji/rji;
%
%  Find j->k vector
%
xki = xk - xj;
yki = yk - yj;
zki = zk - zj;

%  Find unit vector normal to the i-j-k plane

xn = yji*zki - yki*zji;
yn = zji*xki - zki*xji;
zn = xji*yki - xki*yji;
rn = xn*xn + yn*yn + zn*zn;
rn = sqrt(rn);
xn = xn/rn;
yn = yn/rn;
zn = zn/rn;
%  Find unit vector normal to the other 2 directions already found

%  Since original vectors are normalised the result should be likewise

xp = yn*zji - yji*zn;
yp = zn*xji - zji*xn;
zp = xn*yji - xji*yn;

%  Find distances along each unit vector

rji = r*cos(theta);
rn  = r*sin(theta)*sin(phi);
rp  = r*sin(theta)*cos(phi);

%  Multiply unit vector by distances and add to origin to get position

x = xi - rji*xji + rn*xn + rp*xp;
y = yi - rji*yji + rn*yn + rp*yp;
z = zi - rji*zji + rn*zn + rp*zp;
%if abs(norm([x, y, z] - [xi, yi, zi]) - r) > 0.05
%    disp('wrong x y z')
%    [x, y ,z]
%    r
%end

