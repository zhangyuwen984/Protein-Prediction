function [new_moles,deviation]= move_along_SoftMode_molMutation(molecules, eigvector, MtypeLIST)

global ORG_STRUC

vec = zeros(1,3);
new_moles = molecules;
translations=zeros(length(MtypeLIST),3);
rotations=zeros(length(MtypeLIST),3);
count=1;
for i=1:length(MtypeLIST)

        numatom=length(ORG_STRUC.STDMOL(MtypeLIST(i)).types);
        total_inertia=0;
    for j=1:numatom

        dis=[eigvector((count-1)*3+1) eigvector((count-1)*3+2) eigvector((count-1)*3+3)];
        new_moles(i).MOLCOORS(j,1)=molecules(i).MOLCOORS(j,1)-molecules(i).MOLCENTER(1);
        new_moles(i).MOLCOORS(j,2)=molecules(i).MOLCOORS(j,2)-molecules(i).MOLCENTER(2);
        new_moles(i).MOLCOORS(j,3)=molecules(i).MOLCOORS(j,3)-molecules(i).MOLCENTER(3);
%%%%displacement: pararrel and vectical 
        total_inertia=total_inertia+(norm(new_moles(i).MOLCOORS(j,:)))^2;
        direction=new_moles(i).MOLCOORS(j,:);
        if (norm(new_moles(i).MOLCOORS(j,:))<0.05)
            translation(j,:)=dis(:)';
            rotation(j,:)=[0 0 0];
        else
            translation(j,:)=dot(dis,direction)*direction/(norm(direction)*norm(direction));
            rotation(j,:)=cross(dis(:)',direction);
        end
        translations(i,:)=translations(i,:)+translation(j,:);
        rotations(i,:)=rotations(i,:)+rotation(j,:);
        count=count+1;
    end
        translations(i,:) = translations(i,:)/numatom;
        if numatom > 1
           rotations(i,:) = rotations(i,:)/total_inertia;
        end
        normrot(i) = norm(rotations(i,:));
end
for i=1:length(MtypeLIST)
     numatom=length(ORG_STRUC.STDMOL(MtypeLIST(i)).types);
     if numatom > 1
         rotations(i,:)=rotations(i,:)/normrot(i);
         x=rotations(i,1);
         y=rotations(i,2);
         z=rotations(i,3);
         theta=normrot(i);             %translation
         M=[cos(theta)+(1-cos(theta))*(x^2) (1-cos(theta))*x*y-sin(theta)*z (1-cos(theta))*x*z+sin(theta)*y; ...
            (1-cos(theta))*y*x+sin(theta)*z cos(theta)+(1-cos(theta))*(y^2) (1-cos(theta))*y*z-sin(theta)*x; ...
            (1-cos(theta))*z*x-sin(theta)*y  (1-cos(theta))*z*y+sin(theta)*x cos(theta)+(1-cos(theta))*(z^2) ];
         new_moles(i).MOLCOORS=new_moles(i).MOLCOORS*M;
     else
         new_moles(i).MOLCOORS=new_moles(i).MOLCOORS;  %no rotations for single atom
     end
            new_moles(i).MOLCENTER = molecules(i).MOLCENTER + translations(i,:);
    for j=1:numatom
        new_moles(i).MOLCOORS(j,:) = new_moles(i).MOLCOORS(j,:)+new_moles(i).MOLCENTER;
        deviation(i,j,:) = new_moles(i).MOLCOORS(j,:) -molecules(i).MOLCOORS(j,:);
    end
    new_moles(i).ZMATRIX = real(NEW_coord2Zmatrix(new_moles(i).MOLCOORS, ORG_STRUC.STDMOL(MtypeLIST(i)).format));
end

